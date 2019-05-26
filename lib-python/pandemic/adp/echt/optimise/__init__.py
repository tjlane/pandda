import os, copy, collections
from libtbx import adopt_init_args, group_args
from bamboo.common.logs import Log
import numpy

def optimisation_summary_from_target_evaluator(target_evaluator):

        te = target_evaluator

        """Print the optimisation values and weights"""
        s = te.log._bar()
        s += '\nOptimisation Summary: {}\n'.format(self.label)
        s += te.log._bar()
        s += '\nOptimisation Least Squares: {:5.3f}'.format(te.optimisation_lst_sqr)
        s += '\nOptimisation Penalty: {:5.3f}'.format(te.optimisation_penalty)
        s += '\nOptimisation Target: {:5.3f}'.format(te.optimisation_lst_sqr+te.optimisation_penalty)
        s += '\n'+self.log._bar()
        s += '\nModels used for optimisation (and total weights):'
        wts = self._op.wgts.mean(axis=1)
        assert len(wts) == len(self.get_dataset_mask())
        for i_rel, i_abs in enumerate(self.get_dataset_mask()):
            s += '\n\t{:>3d} -> {:.3f}'.format(i_abs,wts[i_rel])
        s += '\n'+te.log._bar()
        if show: te.log(s)
        return s


class LevelTargetUijCalculator:


    def __init__(self,
            uij_target,
            ):
        adopt_init_args(self, locals())

    def __call__(self,
            uij_fitted,
            ignore_level=None,
            ):
        arr = uij_fitted.copy()
        # Zero the current level (as should not be included in target function!)
        if ignore_level is not None:
            arr[ignore_level] = 0.0
        # Sum over all other levels to find total currently modelled
        arr_sum = numpy.sum(arr, axis=0)
        assert arr_sum.shape == self.uij_target.shape
        # Subtract currently modelled from input data
        return self.uij_target - arr_sum


class OptimiseEchtModel:


    def __init__(self,
            optimise_tls_function,
            optimise_adp_function,
            optimise_level_amplitudes_function,
            n_cpus = 1,
            verbose = False,
            log = None,
            ):
        if log is None: log = Log()
        adopt_init_args(self, locals())

        from pandemic.adp.echt.optimise.tls import OptimiseTLSLevel
        self.optimise_tls_level = OptimiseTLSLevel(
            optimisation_function = self.optimise_tls_function,
            n_cpus = n_cpus,
            )

        if self.optimise_adp_function is not None:
            from pandemic.adp.echt.optimise.uij import OptimiseUijLevel
            self.optimise_adp_level = OptimiseUijLevel(
                optimisation_function = self.optimise_adp_function,
                n_cpus = n_cpus,
                )
        else:
            self.optimise_adp_level = None

        self.optimise_level_amplitudes = optimise_level_amplitudes_function

        from pandemic.adp.echt.optimise.repair import RepairEchtModel
        self.repair_model = RepairEchtModel()

    def __call__(self,
            model_object,
            level_group_connections,
            uij_target,
            uij_target_weights,
            uij_isotropic_mask,
            tracking_object,
            uij_optimisation_mask = None,
            ):

        self.log('\nOptimising TLS models for {} levels'.format(model_object.n_tls_levels) \
            + ' (+ {} level)'.format(model_object.adp_level_name)*(self.optimise_adp_level is not None))

        # Increment counter at beginning of cycle (ensure we have a new number)
        tracking_object.i_cycle += 1

        # Class for calculating target for each optimisation
        calculate_uij_target = LevelTargetUijCalculator(uij_target)

        tracking_object.update(
            uij_target = uij_target,
            uij_lvl = model_object.uijs(),
            step = 'start',
            )

        # Iterate through the TLS levels of the fitting
        for i_level in xrange(model_object.n_tls_levels):

            self.log.subheading('Macrocycle {}: '.format(tracking_object.i_cycle)+'Fitting TLS Groups (level {} - {})'.format(i_level+1, model_object.tls_level_names[i_level]))

            # Optimise the groups for one level
            model_object.tls_objects[i_level] = self.optimise_tls_level(
                uij_target = calculate_uij_target(uij_fitted=model_object.uijs(), ignore_level=i_level),
                uij_target_weights = uij_target_weights,
                uij_isotropic_mask = uij_isotropic_mask,
                uij_optimisation_mask = uij_optimisation_mask,
                tls_groups = model_object.tls_objects[i_level],
                tls_selections = model_object.tls_selections[i_level],
                )

            # Optimise the amplitudes between levels
            if self.optimise_level_amplitudes is not None:

                self.log.subheading('Macrocycle {}: '.format(tracking_object.i_cycle)+'Optimising inter-level amplitudes')

                # Ensure no negative eigenvalues
                self.repair_model(model_object)

                for max_recursions in [1, None]:
                    if max_recursions is None:
                        self.log('Optimising all levels simultaneously')
                    else:
                        self.log('Optimising sets of {} neighbouring levels'.format(max_recursions+1))
                    self.optimise_level_amplitudes(
                        uij_target = uij_target,
                        uij_target_weights = uij_target_weights,
                        uij_isotropic_mask = uij_isotropic_mask,
                        model_object = model_object,
                        level_group_connections  = level_group_connections,
                        max_recursions = max_recursions,
                        # dataset mask
                        )

            # Update tracking
            tracking_object.update(
                uij_target = uij_target,
                uij_lvl = model_object.uijs(),
                step = 'level {}'.format(i_level+1),
                i_level = i_level,
                write_graphs = False,
                )

        if self.optimise_adp_level is not None:

            # Fit the residuals
            self.log.subheading('Macrocycle {}: '.format(tracking_object.i_cycle)+'Optimising residual atomic Uijs')

            # Update the target uij by subtracting contributions from other levels
            model_object.adp_values = self.optimise_adp_level(
                uij_values = model_object.adp_values,
                uij_target = calculate_uij_target(uij_fitted=model_object.uijs(), ignore_level=-1),
                uij_target_weights = uij_target_weights,
                uij_isotropic_mask = uij_isotropic_mask,
                uij_optimisation_mask = uij_optimisation_mask,
                )

            # Optimise the amplitudes between levels
            if self.optimise_level_amplitudes is not None:

                self.log.subheading('Macrocycle {}: '.format(tracking_object.i_cycle)+'Optimising inter-level amplitudes')

                # Ensure no negative eigenvalues
                self.repair_model(model_object)

                for max_recursions in [1, None]:
                    if max_recursions is None:
                        self.log('Optimising all levels simultaneously')
                    else:
                        self.log('Optimising sets of {} neighbouring levels'.format(max_recursions+1))
                    self.optimise_level_amplitudes(
                        uij_target = uij_target,
                        uij_target_weights = uij_target_weights,
                        uij_isotropic_mask = uij_isotropic_mask,
                        model_object = model_object,
                        level_group_connections  = level_group_connections,
                        max_recursions = max_recursions,
                        )

            # Update tracking
            tracking_object.update(
                uij_target = uij_target,
                uij_lvl = model_object.uijs(),
                step = 'residual',
                i_level = model_object.n_levels-1,
                write_graphs = False,
                )

        # Ensure no negative eigenvalues
        self.repair_model(model_object)

        # Update tracking
        tracking_object.update(
            uij_target = uij_target,
            uij_lvl = model_object.uijs(),
            step = 'inter-level',
            i_level = range(model_object.n_levels),
            write_graphs = True,
            )

        return model_object


class UpdateOptimisationFunction:


    level_amplitude_string = 'level amplitudes weights'

    def __init__(self,
        n_cycles,
        gradient_optimisation_decay_factor,
        model_optimisation_function,
        plotting_object = None,
        verbose = True,
        log = None,
        ):
        if log is None: log = Log()
        assert n_cycles >= 1
        start_values = group_args(
            level_amplitudes_optimisation_weights = copy.deepcopy(model_optimisation_function.optimise_level_amplitudes.optimisation_weights),
            )
        history = collections.OrderedDict()
        adopt_init_args(self, locals(), exclude=('model_optimisation_function',))

    def update(self,
        model_optimisation_function,
        i_cycle,
        ):

        assert i_cycle >= 0
        assert i_cycle < self.n_cycles

        # Percentage of program remaining
        if self.n_cycles <= 1:
            frac = 0.0
        else:
            frac = float(self.n_cycles-i_cycle) / float(self.n_cycles)

        decay_factor = self.gradient_optimisation_decay_factor
        total_decay_factor = decay_factor ** i_cycle

        # Update parameters
        self.log.subheading('Updating Level Amplitude Optimisation Weights')
        self.log('> Cycle {}'.format(i_cycle+1))
        self.log('> Total decay factor (relative to starting values): {}'.format(total_decay_factor))
        for k, v in self.start_values.level_amplitudes_optimisation_weights.__dict__.iteritems():
            new_v = v / total_decay_factor
            self.log('Updating {} = {} -> {}'.format(k, v, new_v))
            model_optimisation_function.optimise_level_amplitudes.optimisation_weights.__dict__[k] = new_v
            self.history \
                .setdefault(self.level_amplitude_string, collections.OrderedDict()) \
                .setdefault(k, []).append((i_cycle+1, new_v))

    def write(self, output_directory):

        output_files = {}

        for s in [self.level_amplitude_string]:
            prefix = os.path.join(output_directory, s.replace(' ', '_'))
            for variable, values in self.history[s].iteritems():
                x_vals, y_vals = zip(*values)
                filename = prefix+'-{variable}.png'.format(variable=variable)
                fig, axis = self.plotting_object.lineplot(
                    x_vals = x_vals,
                    y_vals = y_vals,
                    title = 'Weight "{}" over cycles'.format(variable),
                    x_label = 'cycle',
                    y_label = 'weight',
                    marker = '*',
                    filename = None, # returns fig and axis
                    )
                axis.set_xticks(x_vals)
                axis.set_xticklabels(map(int,x_vals))
                axis.set_yscale('log', base=10.)
                self.plotting_object.helper.write_and_close_fig(fig=fig, filename=filename)
                output_files.setdefault(s,collections.OrderedDict())[variable] = filename

        return output_files
