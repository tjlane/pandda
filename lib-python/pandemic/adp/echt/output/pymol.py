import giant.logs as lg
logger = lg.getLogger(__name__)

import os, collections

from giant.structure.pymol import auto_chain_images


class WriteEchtPymolImages(object):

    output_key = 'level_uijs_pymol_by_chain_png'
    output_path_template_prefix = 'level_{level_num}'

    def __init__(self, output_directory):
        
        self.output_directory = output_directory

        assert '{level_num}' in self.output_path_template_prefix

    def __call__(self,
        level_pdbs_dict,
        model_object,
        reference_hierarchy,
        overall_atom_mask,
        ):

        logger.subheading('Generating pymol images for levels in ECHT model')

        if not self.output_directory.exists():
            self.output_directory.mkdir(parents=True)

        ###

        output_files = collections.OrderedDict()

        ###

        # Record any failed/missing files
        missing_files = []

        for i_level, level_name in enumerate(model_object.all_level_names):

            f_prefix = str(
                self.output_directory / self.output_path_template_prefix.format(
                    level_num = i_level+1,
                    )
                )

            of = auto_chain_images(
                structure_filename = level_pdbs_dict[level_name],
                output_prefix = f_prefix,
                style = 'lines+ellipsoids',
                colours = 'bfactor',
                width=1000, height=750,
                )

            if (not of):
                logger.warning('no images have been generated: {}...'.format(f_prefix))
            else:
                missing_files.extend([v for v in of.values() if not os.path.exists(v)])
                output_files[level_name] = of

        if missing_files:
            logger.warning(
                '\n'.join(
                    ['image does not exist: {}'.format(v) for v in missing_files]
                    )
                )

        return {self.output_key : output_files}


class WriteEchtPymolScript(object):

    output_key = 'pymol_script'
    output_path = 'pymol_script.py'

    def __init__(self, output_directory):
        
        self.output_directory = output_directory

    def __call__(self,
        file_dict,
        ):

        from giant.pymol_utils import PythonScript

        if not self.output_directory.exists():
            self.output_directory.mkdir(parents=True)

        s = PythonScript(
            pretty_models = False, 
            pretty_maps = False,
            )

        s.change_into_directory_maybe(
            path = os.path.abspath(str(self.output_directory)),
            )

        for f in [file_dict.get('target_uijs_pdb',None)] + list(file_dict.get('level_uijs_pdb',{}).values()):
            if f is None: continue
            obj = os.path.basename(f)
            s.load_pdb(
                f_name = os.path.relpath(os.path.abspath(f), start=str(self.output_directory)),
                obj = obj,
                )
            s.custom('spectrum', expression='b', selection=obj)
            
        for f in [file_dict.get('output_uijs_pdb',None)] + list(file_dict.get('level_uijs_by_mode_pdb',{}).values()):
            if f is None: continue
            obj = os.path.basename(f)
            s.load_pdb(
                f_name = os.path.relpath(os.path.abspath(f), start=str(self.output_directory)),
                obj = obj,
                )
            s.disable(
                obj = obj,
                )
            s.custom('spectrum', expression='b', selection=obj)

        s.set('ellipsoid_probability', '0.95')
        s.show_as(obj='all', style='lines')
        s.show(obj='all', style='ellipsoids')
        # s.custom('spectrum', expression='b', selection='all')
        s.set('grid_mode', 1)

        s.orient(obj='all')

        filename = str(
            self.output_directory / self.output_path
            )

        s.write_script(filename)

        return {self.output_key : filename}


