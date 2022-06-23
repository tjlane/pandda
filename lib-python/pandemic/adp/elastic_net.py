import giant.logs as lg
logger = lg.getLogger(__name__)

import os, copy, collections
from libtbx import adopt_init_args, group_args

from mmtbx.tls.optimise_amplitudes import OptimisationWeights


class ReplaceableOptimisationWeights(OptimisationWeights):
    """
    Extended weights object for OptimiseAmplitudes class.

    Parameters
    ----------
    sum_of_amplitudes : float
      weight for the lasso-like term (constrains sum of amplitudes)
    sum_of_squared_amplitudes : float
      weight for the ridge-regression-like term (restrains each amplitude to zero)
    sum_of_amplitudes_squared : float
      weight for the square of the lasso-like term (restrains sum of amplitudes to zero)
    """

    @classmethod
    def defaults(cls):
        """Initialise class with default values for each of the fields"""
        return cls(**{f:0.0 for f in cls._fields})

    def transfer_from_other(self, other, require_all=False):
        """
        Transfer weights from an input object (with equivalent attributes or indexed)
        Returns a new object using _replace function of namedtuple
        """

        update_dict = {}

        for f in self._fields:

            # Extract value for the field name (if available or required)
            try:
                v = getattr(other, f)
            except AttributeError as e:
                if require_all is True:
                    raise ValueError('{}\n`other` does not have attribute or contain item: {}'.format(str(e), f))
                else:
                    continue

            # store temporarily in dictionary
            update_dict[f] = v

        if len(update_dict) == 0:
            raise ValueError('No overlapping attributes/items were found between self and `other`')

        return self._replace(**update_dict)


class UpdateOptimisationFunction(object):

    level_amplitude_string = 'level amplitudes weights'

    def __init__(self,
        initial_weights,
        weight_decay_factor,
        weights_to_update,
        minimum_weight,
        output_directory,
        plotting_object = None,
        ):
        # Copy to be ensure detached object
        initial_weights = copy.deepcopy(initial_weights)
        history = collections.OrderedDict()
        result = group_args(
            output_files = None,
            )
        adopt_init_args(self, locals())

    def get_weights(self, model_optimisation_function):
        return model_optimisation_function.optimise_level_amplitudes.optimisation_weights

    def set_weights(self, model_optimisation_function, weights):
        model_optimisation_function.optimise_level_amplitudes.optimisation_weights = weights

    def update(self,
        model_optimisation_function,
        n_cycle,
        ):

        assert n_cycle >= 1

        start_weights = self.initial_weights
        minimum_weight = self.minimum_weight
        decay_factor = self.weight_decay_factor
        total_decay_factor = decay_factor ** (n_cycle-1)

        # Update parameters
        logger.subheading('Updating Level Amplitude Optimisation Weights')
        logger('> Cycle {}\n'.format(n_cycle))
        logger('> Total decay factor (relative to starting values): {}\n'.format(total_decay_factor))

        new_weights_dict = {}
        for k in self.weights_to_update:
            sta_v = getattr(start_weights, k)
            new_v = sta_v * total_decay_factor
            logger('Updating {} = {} -> {}'.format(k, sta_v, new_v))
            new_weights_dict[k] = new_v
        # Update weights if above threshold
        elastic_net_weight_total = (
            new_weights_dict['sum_of_amplitudes'] + 
            new_weights_dict['sum_of_squared_amplitudes']
            )
        if (minimum_weight is not None) and (elastic_net_weight_total < minimum_weight):
            logger(
                (
                    '\n*** Sum of elastic-net optimisation weights ({}) '
                    'is below the minimum threshold ({}). Not updating weights. ***'
                    ).format(elastic_net_weight_total, minimum_weight)
                )
        else:
            # Get current weight object
            current_weights = self.get_weights(model_optimisation_function)
            # Create new object, overriding the selected weights
            new_weights = current_weights.transfer_from_other(group_args(**new_weights_dict))
            # Replace the original weights object with the new weights object
            self.set_weights(model_optimisation_function, new_weights)

        # (re)extract current weights
        current_weights = self.get_weights(model_optimisation_function)

        # Update history
        logger('\nCurrent values:')
        for k, v in current_weights._asdict().items():
            logger('{} -> {}'.format(k, v))
            self.history \
                .setdefault(self.level_amplitude_string, collections.OrderedDict()) \
                .setdefault(k, []).append((n_cycle, v))

    def write_output(self):

        logger.subheading('Writing optimisation weights graphs')

        output_files = collections.OrderedDict()

        for s in [self.level_amplitude_string]:
            prefix = os.path.join(self.output_directory, s.replace(' ', '_'))
            for variable, values in self.history[s].items():
                x_vals, y_vals = list(zip(*values))
                filename = prefix+'-{variable}.png'.format(variable=variable)
                fig, axis = self.plotting_object.lineplot(
                    x_vals = x_vals,
                    y_vals = y_vals,
                    title = 'Weight "{}" over cycles'.format(variable),
                    x_label = 'cycle',
                    y_label = 'weight',
                    x_ticks = list(map(int,x_vals)),
                    filename = None, # returns fig and axis
                    background_line_type = 'chunky',
                    )
                # Log scale if data has positive values
                if max(y_vals) > 0.0:
                    try:
                        axis.set_yscale('log')
                    except:
                        pass
                self.plotting_object.helper.write_and_close_fig(fig=fig, filename=filename)
                output_files.setdefault(s,collections.OrderedDict())[variable] = filename

        self.result.output_files = output_files

        return output_files

    def as_html_summary(self):
        from pandemic.adp.echt.html.optimise import EchtOptimisationParametersHtmlSummary
        return EchtOptimisationParametersHtmlSummary(self)
