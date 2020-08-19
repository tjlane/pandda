import giant.logs as lg
logger = lg.getLogger(__name__)

import os, collections
from libtbx import adopt_init_args
from giant.paths import easy_directory

from scitbx.array_family import flex

from giant.structure.uij import uij_to_b


class WriteEchtStructures:


    pymol_script_py = 'pymol_script.py'

    def __init__(self,
        output_directory,
        ):
        adopt_init_args(self, locals())

    def show(self, label, pdbs, max=5):
        logger(label)
        for l, p in list(pdbs.iteritems())[:max]+[('', '...')]*(len(pdbs)>max):
            logger('> {}: {}'.format(l, p))

    def __call__(self,
        level_group_array,
        model_object,
        isotropic_mask,
        models,
        overall_selection = None,
        overall_selection_bool = None,
        ):
        """Extract output and generate summaries"""

        from pandemic.adp.echt.output.tls_headers import MultiModelTLSHeaderFactory
        make_tls_headers = MultiModelTLSHeaderFactory(
            tls_group_array = level_group_array,
            tls_selection_strings = model_object.tls_selection_strings,
            tls_objects = model_object.tls_objects,
            overall_selection = overall_selection,
            )

        from pandemic.adp.output import MultiModelStructureWriter
        output_structures = MultiModelStructureWriter(
            output_directory = self.output_directory,
            models = models,
            atom_mask = overall_selection_bool,
            isotropic_mask = isotropic_mask,
            )

        #------------------------------------------------------------------------------#
        #---#                            Extract uijs                              #---#
        #------------------------------------------------------------------------------#

        uij_lvl = model_object.uijs()
        uij_all = uij_lvl.sum(axis=0)
        uij_tls = uij_lvl[:-1].sum(axis=0)

        #------------------------------------------------------------------------------#
        #---#                     Write output structures                          #---#
        #------------------------------------------------------------------------------#

        output_dict = {}

        # All TLS levels
        tls_headers = make_tls_headers(i_levels=None)

        # Fully parameterised structures
        logger.subheading('Writing full uij models')
        # Full models (with original values for non-fitted atoms)
        pdbs = output_structures(
            uij = uij_all,
            iso = map(uij_to_b, uij_all),
            headers = tls_headers,
            model_suffix = '.all.pdb',
            blank_copy = False)
        self.show('All levels integrated into original structures', pdbs)
        output_dict['complete_structures'] = pdbs

        # Full models (zero values for non-fitted atoms)
        pdbs = output_structures(
            uij = uij_all,
            iso = map(uij_to_b, uij_all),
            headers = tls_headers,
            model_suffix = '.all-levels.pdb',
            blank_copy = True)
        self.show('All levels', pdbs)
        output_dict['all_components'] = pdbs

        # TLS-parameterised structures (total TLS)
        logger.subheading('Writing all TLS contributions')
        pdbs = output_structures(
            uij = uij_tls,
            iso = map(uij_to_b, uij_tls),
            headers = tls_headers,
            model_suffix = '.all-tls-levels.pdb',
            blank_copy = True)
        self.show('All tls levels', pdbs)
        output_dict['all_tls_levels'] = pdbs

        # Atomic level contributions
        logger.subheading('Writing atomic level')
        uij_atom = uij_lvl[-1]
        pdbs = output_structures(
            uij = uij_atom,
            iso = map(uij_to_b, uij_atom),
            headers = None,
            model_suffix = '.atomic-level.pdb',
            blank_copy = True)
        self.show('{} level'.format(model_object.adp_level_name.title()), pdbs)
        output_dict['atomic_level'] = pdbs

        # Level by level TLS-parameterised structures (single level contribution)
        logger.subheading('Writing individual TLS levels')
        for i_level in xrange(model_object.n_tls_levels):
            uij_this = uij_lvl[i_level]
            pdbs = output_structures(
                uij = uij_this,
                iso = map(uij_to_b, uij_this),
                headers = make_tls_headers(i_levels=[i_level]),
                model_suffix = '.tls-level-{:04}.pdb'.format(i_level+1),
                blank_copy = True)
            self.show('Level {}'.format(i_level+1), pdbs)
            output_dict.setdefault('tls_levels', collections.OrderedDict())[i_level] = pdbs

        # Level by level TLS-parameterised structures (cumulative level contributions)
        logger.subheading('Writing different combinations of TLS levels')
        for i_level, j_level in flex.nested_loop((model_object.n_tls_levels, model_object.n_tls_levels)):
            if i_level >= j_level: continue
            cuml_uij = uij_lvl[i_level:j_level+1].sum(axis=0)
            pdbs = output_structures(
                uij = cuml_uij,
                iso = map(uij_to_b, cuml_uij),
                headers = make_tls_headers(i_levels=range(i_level, j_level+1)),
                model_suffix = '.tls-level-{:04}-to-{:04}.pdb'.format(i_level+1,j_level+1),
                blank_copy = True)
            self.show('Levels {}-{}'.format(i_level+1, j_level+1), pdbs)
            output_dict.setdefault('tls_levels', collections.OrderedDict())[(i_level,j_level)] = pdbs

        logger.subheading('Writing convenience pymol scripts')
        pymol_scripts = self.write_pymol_scripts(
            model_keys = [m.tag for m in models],
            file_dict = output_dict,
            )
        self.show('Pymol scripts for viewing each structure', pymol_scripts)
        output_dict['pymol_scripts'] = pymol_scripts

        return output_dict

    def write_pymol_scripts(self, model_keys, file_dict):

        from giant.pymol_utils import PythonScript

        output_files = collections.OrderedDict()

        for key in model_keys:

            # Output directory and filename
            py_d = easy_directory(os.path.join(self.output_directory, key))
            py_f = os.path.join(py_d, self.pymol_script_py)

            s = PythonScript(pretty_models=False, pretty_maps=False)

            s.change_into_directory_maybe(path=os.path.abspath(py_d))

            f_list = []

            f = file_dict.get('all_components',{}).get(key)
            if (f is not None): f_list.append(f)

            for i_l, f_dict in file_dict.get('tls_levels',{}).iteritems():
                f = f_dict.get(key)
                if (f is not None):
                    f_list.append(f)

            f = file_dict.get('atomic_level',{}).get(key)
            if (f is not None): f_list.append(f)

            for f in f_list:
                obj = os.path.basename(f)
                s.load_pdb(f_name=os.path.relpath(os.path.abspath(f), start=py_d), obj=obj)

            s.show_as(obj='all', style='lines')
            s.show(obj='all', style='ellipsoids')
            s.custom('spectrum', expression='b', selection='all')
            s.set('grid_mode', 1)

            s.orient(obj='all')

            s.write_script(py_f)

            output_files[key] = py_f

        return output_files
