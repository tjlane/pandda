import giant.logs as lg
logger = lg.getLogger(__name__)

import os, collections
from libtbx import adopt_init_args
from giant.paths import easy_directory

import numpy
from scitbx.array_family import flex

from pandemic.adp.hierarchy.utils import StructureFactory


class MultiModelStructureWriter(object):


    def __init__(self,
        output_directory,
        models,
        atom_mask,
        isotropic_mask,
        ):

        # Make sure output directory exists
        output_directory = easy_directory(output_directory)

        atom_mask = flex.bool(atom_mask)

        structure_factories = [StructureFactory(master_h=m.hierarchy) for m in models]

        # Extract useful integers
        n_datasets = len(models)
        n_atoms = sum(atom_mask)

        adopt_init_args(self, locals())

    def __call__(self,
        uij,
        iso = None,
        headers = None,
        model_suffix = '.pdb',
        blank_copy = False,
        ):

        n_datasets = self.n_datasets
        n_atoms = self.n_atoms

        # Validate AnoUijs
        uij = numpy.array(uij)
        assert uij.shape == (n_datasets, n_atoms, 6)

        # Validate IsoBs
        if iso is not None:
            iso = numpy.array(iso)
            assert iso.shape == (n_datasets, n_atoms)

        # Validate header lines
        if headers is not None:
            if isinstance(headers[0], str):
                headers = [headers]
            assert set(map(len,headers)) == {len(self.models)}
            open_append = True
        else:
            open_append = False

        # List of output pdbs
        pdbs = collections.OrderedDict()

        # Apply to each model and write out
        for i, mdl in enumerate(self.models):

            # Create model paths
            mdl_d = easy_directory(os.path.join(self.output_directory, mdl.tag))
            mdl_f = os.path.join(mdl_d, mdl.tag+model_suffix)

            # Store in output dict
            pdbs[mdl.tag] = mdl_f

            try:
                # Create new copy of hierarchy and extract atoms
                h = self.structure_factories[i].custom_copy(
                    uij = self.isotropic_mask(uij[i]),
                    iso = (iso[i] if (iso is not None) else None),
                    mask = self.atom_mask,
                    blank_copy = blank_copy,
                    )

                # Write headers
                if headers is not None:
                    with open(mdl_f, 'w') as fh:
                        [fh.write(l[i].strip()+'\n') for l in headers]

                # Write atoms
                crystal_symmetry = (mdl.crystal_symmetry if hasattr(mdl, 'crystal_symmetry') else None)
                h.write_pdb_file(
                        mdl_f,
                        open_append = open_append,
                        crystal_symmetry = crystal_symmetry,
                        )

            except Exception as e:
                import traceback
                msg = 'Failed to write structure for dataset {}: ({})\n\n{}'.format(mdl.tag, mdl_f, traceback.format_exc())
                logger.warning(msg)
                continue

        return pdbs


