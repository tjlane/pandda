import giant.logs as lg
logger = lg.getLogger(__name__)

from libtbx import adopt_init_args, group_args
from giant.exceptions import Sorry, Failure
import numpy

def scale_weights(weights, scale):

    import copy
    new_weights = copy.deepcopy(weights)

    for k in new_weights.__dict__:
        if k.startswith('_'): continue
        new_weights.__dict__[k] = scale * new_weights.__dict__[k]

    return new_weights

def combine_optimisation_weights(
        atom_weight_array,
        dataset_weight_array,
        ):
    # Calculate combined weighting
    total_weight_array = dataset_weight_array.reshape((len(dataset_weight_array),1)) * atom_weight_array
    # Normalise the total weights
    total_weight_array /= total_weight_array.mean()
    assert abs(total_weight_array.mean()-1.0) < 1e-3, 'average weight should be 1!'
    return total_weight_array


class AtomWeightCalculator(object):


    _power_hash = {
            'one' : 0.0,
            'inverse_mod_U' : -1.0,
            'inverse_mod_U_squared' : -2.0,
            'inverse_mod_U_cubed' : -3.0,
            }

    def __init__(self,
            weighting,
            renormalise_by_dataset = True,
            ):
        assert weighting in self._power_hash
        power = self._power_hash[weighting]
        adopt_init_args(self, locals())

    def __call__(self,
            uij_array,
            dataset_labels = None,
            ):

        if dataset_labels is None:
            dataset_labels = list(range(1, len(uij_array)+1))

        # Calculate isotropic equivalent
        uij_modulus = uij_array[:,:,0:3].mean(axis=2)

        # Check to see if any u have zero size
        if (uij_modulus == 0.0).any():
            zero_b_atoms = list(zip(*numpy.where((uij_modulus == 0.0))))
            message = 'Some atoms have zero-value b-factors!'
            message += '\n\t'+'\n\t'.join(['atom {} in dataset {}'.format(i_a,dataset_labels[i_d]) for i_d,i_a in zero_b_atoms])
            if (self.weighting == 'one'):
                logger.warning(message)
            else:
                raise Sorry('Cannot calculate atom weights: '+message)

        # Calculate weights
        weights = numpy.power(uij_modulus, self.power)

        # Normalise individually by dataset
        if self.renormalise_by_dataset:
            for i,c in enumerate(weights):
                weights[i] = c / c.mean()

        # Normalise to average of one
        weights = weights / weights.mean()

        assert abs(weights.mean()-1.0) < 1e-3, 'average weight should be approximately 1!'

        return weights


class DatasetWeightCalculator(object):


    _power_hash = {
            'one' : 0.0,
            'inverse_resolution' : -1.0,
            'inverse_resolution_squared' : -2.0,
            'inverse_resolution_cubed' : -3.0,
            }

    def __init__(self,
            weighting,
            ):
        assert weighting in self._power_hash
        power = self._power_hash[weighting]
        adopt_init_args(self, locals())

    def __call__(self,
            resolutions = None,
            dataset_labels = None,
            ):

        if (resolutions is None) and (dataset_labels is None):
            raise Exception('must provide resolutions or dataset_labels')
        elif (resolutions is None):
            assert self.weighting == 'one'
            resolutions = numpy.ones_like(dataset_labels, dtype=float)
        elif (dataset_labels is None):
            dataset_labels = list(range(1, len(resolutions)+1))

        resolutions = numpy.array(resolutions)
        if not (resolutions > 0.0).all():
            raise Exception('Negative resolutions have been provided: \n{}'.format(resolutions))

        # Create weights
        weights = numpy.power(resolutions, self.power)

        # Normalise
        weights = weights / weights.mean()

        assert abs(weights.mean()-1.0) < 1e-3, 'average weight should be approximately 1!'

        return weights


class UijArrayWeightsTask(object):

    def __init__(self,
            dataset_weighting = 'one',
            atom_weighting = 'one',
            renormalise_by_dataset = False,
            ):

        calculate_dataset_weights = DatasetWeightCalculator(
            weighting = dataset_weighting,
            )

        calculate_atom_weights = AtomWeightCalculator(
            weighting = atom_weighting,
            renormalise_by_dataset = renormalise_by_dataset,
            )

        adopt_init_args(self, locals())

    def run(self,
            uij_values,
            resolutions = None,
            dataset_labels = None,
            ):

        if dataset_labels is None:
            dataset_labels = list(range(1, len(uij_values)+1))

        assert len(dataset_labels) == len(uij_values)

        if resolutions is not None:
            assert len(resolutions) == len(uij_values)

            if resolutions.count(None):
                assert resolutions.count(None) == len(resolutions), 'either all resolutions must be provided or none: {}'.format(resolutions)
                resolutions = None

        dataset_weight_array = self.calculate_dataset_weights(
            resolutions = resolutions,
            dataset_labels = dataset_labels,
            )

        atom_weight_array = self.calculate_atom_weights(
            uij_array = uij_values,
            dataset_labels = dataset_labels,
            )

        total_weight_array = combine_optimisation_weights(
            atom_weight_array = atom_weight_array,
            dataset_weight_array = dataset_weight_array,
            )

        self.result = group_args(
            atom_weight_array = atom_weight_array,
            dataset_weight_array = dataset_weight_array,
            total_weight_array = total_weight_array,
            )

        return self.result


