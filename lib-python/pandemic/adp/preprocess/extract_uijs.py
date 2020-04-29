import giant.logs as lg
logger = lg.getLogger(__name__)

from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure

import numpy
from pandemic.adp import constants


class ExtractAndProcessModelUijsTask:


    def __init__(self,
            expected_disorder_model,
            verbose = False,
            ):
        actual_disorder_model = expected_disorder_model
        adopt_init_args(self, locals())

    @staticmethod
    def extract_uijs_from_models(
            models_atoms,
            overall_atom_mask = None,
            ):
        # Extract all atoms from all datasets
        assert len(models_atoms) > 0, 'No models have been used?!'
        assert len(models_atoms[0]) > 0, 'No atoms have been extracted from models?'
        # Extract uij and xyz from all datasets (and only select atoms we're interested in)
        model_uij = numpy.array([a.extract_uij() for a in models_atoms])
        # Filter on input mask
        if overall_atom_mask is not None:
            model_uij = model_uij[:,overall_atom_mask]
        return model_uij

    def determine_disorder_model_from_mask(self,
            isotropic_mask,
            ):
        # Initialise disorder model as expected
        actual_disorder_model = self.expected_disorder_model
        # Check that the isotropic mask is as expected
        if isotropic_mask.selection.all_eq(True):
            if (self.expected_disorder_model == 'anisotropic'):
                raise Sorry('All atoms selected for fitting have isotropic ADPs and input.input_uij_model is currently set to "{}".\nOptions:'.format(self.expected_disorder_model) + \
                        '\n   1) Re-refine structures with anisotropic ADPs.' + \
                        '\n   2) Set input.input_uij_model=isotropic to allow use of isotropic ADPs.')
            if (self.expected_disorder_model == 'mixed'):
                logger('Input disorder model is set to "{}", but all input atoms have isotropic Uij.'.format(self.expected_disorder_model))
                actual_disorder_model = 'isotropic'
                logger('Updated disorder_model to "{}"'.format(actual_disorder_model))
        elif isotropic_mask.selection.all_eq(False):
            if (self.expected_disorder_model == 'isotropic'):
                raise Sorry('input.input_uij_model is set to "{}" but all atoms selected have anisotropic Uij.'.format(self.expected_disorder_model))
            if (self.expected_disorder_model == 'mixed'):
                logger('Input disorder model is set to "{}", but all input atoms have anisotropic Uij.'.format(self.expected_disorder_model))
                actual_disorder_model = 'anisotropic'
                logger('Updated disorder_model to "{}"'.format(actual_disorder_model))
        else:
            if (self.expected_disorder_model != 'mixed'):
                raise Sorry('Some atoms for fitting (but not all) have anisotropically-refined ADPs and input.input_uij_model is currently set to "{}".\nOptions:'.format(self.expected_disorder_model) + \
                        '\n   1) Re-refine structures so that all atoms selected have {} ADPs and set input.input_uij_model={}.'.format(self.expected_disorder_model, self.expected_disorder_model) + \
                        '\n   2) Change levels.overall_selection so that atoms with anisotropic/isotropic ADPs are not selected.' + \
                        '\n   3) Change input.input_uij_model to "mixed".')

        return actual_disorder_model

    @staticmethod
    def generate_uij_for_isotropic_atoms(
            models_atoms,
            model_uij,
            isotropic_mask,
            overall_atom_mask = None,
            ):

        # Create copy
        import copy
        model_uij = copy.deepcopy(model_uij)

        # Find atoms that are isotropic
        assert (model_uij[:,isotropic_mask.selection,:] == -1).all()
        # Set atoms to zeros
        model_uij[:,isotropic_mask.selection,:] = 0.0
        # Extract B-values and convert to Uijs
        model_b = numpy.array([a.extract_b()/constants.EIGHTPISQ for a in models_atoms])
        # Apply overall mask
        if overall_atom_mask is not None:
            model_b = model_b[:,overall_atom_mask]
        # Extract values for isotropic atoms
        model_b_iso = model_b[:,isotropic_mask.selection]
        # Set Uij values of isotropic atoms
        model_uij[:,isotropic_mask.selection,0] = model_b_iso
        model_uij[:,isotropic_mask.selection,1] = model_b_iso
        model_uij[:,isotropic_mask.selection,2] = model_b_iso

        return model_uij

    def run(self,
            models,
            overall_atom_mask = None,
            ):

        logger.subheading('Extracting Uijs from models')

        # Extract atoms from models
        models_atoms = [m.hierarchy.atoms() for m in models]

        model_uij = self.extract_uijs_from_models(
                models_atoms = models_atoms,
                overall_atom_mask = overall_atom_mask,
                )

        # Create a mask object for isotropic atoms
        from pandemic.adp.uijs.masks import UijIsotropicMask
        isotropic_mask = UijIsotropicMask.from_uij_array(model_uij, axis=(0,2))

        # Determine the disorder model
        actual_disorder_model = self.determine_disorder_model_from_mask(isotropic_mask)

        # Extract b-values for isotropic atoms and convert to isotropic Uij
        if not isotropic_mask.all_anisotropic:
            model_uij = self.generate_uij_for_isotropic_atoms(
                    models_atoms = models_atoms,
                    model_uij = model_uij,
                    isotropic_mask = isotropic_mask,
                    overall_atom_mask = overall_atom_mask,
                    )

        # All -1 values should now have been removed from input array
        assert not (model_uij==-1).all(axis=2).any()

        self.result = group_args(
                disorder_model = actual_disorder_model,
                model_uij = model_uij,
                isotropic_mask = isotropic_mask,
                )

        self.show_summary()

        return self.result

    def show_summary(self):
        r = self.result

        logger.bar()
        logger('Expected disorder model was {}'.format(self.expected_disorder_model))
        logger('Actual disorder model is {}'.format(r.disorder_model))
        logger.bar()
        logger('Anisotropic atoms: {}'.format(r.isotropic_mask.n_anisotropic))
        logger('Isotropic atoms: {}'.format(r.isotropic_mask.n_isotropic))
        logger.bar()

