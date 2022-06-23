import giant.logs as lg
logger = lg.getLogger(__name__)

import os, glob, collections
from libtbx import adopt_init_args, group_args
from giant.structure.pymol import auto_chain_images
from pandemic.adp.utils import show_file_dict
import numpy

def translate_phenix_selections_to_pymol_selections_simple(selections):
    """Convert simplex phenix selections into simplex pymol selections"""

    easy_regexes = [
            ("chain '?[A-Za-z]{1,2}'?",
                lambda x: x
                ),
            ("chain '?[A-Za-z]{1,2}'? and resid '?[0-9]*?[A-Z]?'? through '?[0-9]*[A-Z]?'?",
                lambda x: x.replace("resid", "resi").replace(" through ",":")
                ),
                    ]

    import re
    output_selections = []
    for s in selections:
        # Start
        o = None
        # Remove double spaces and padding
        while '  ' in s:
            s = s.replace('  ',' ')
        s = s.strip(' ')
        logger.debug('Phenix Selection String: {}'.format(s))
        for rgx_s, trans_func in easy_regexes:
            rgx = re.compile(rgx_s)
            mat = rgx.findall(s)
            # No matches or too many matches
            if len(mat) != 1:
                logger.debug('> No matches or multiple matches to {}'.format(rgx_s))
                continue
            # If the match is the same as input string, then process
            if mat[0] == s:
                o = trans_func(s)
                logger.debug('Exact match for {}'.format(rgx_s))
                logger.debug('Translating to pymol: \n\t{} -> {}'.format(s, o))
                break
        # Append processed selection or none
        output_selections.append(o)

    assert len(output_selections) == len(selections)
    return output_selections


class HierarchySummaryOutput(object):

    def __init__(self,
        output_files,
        ):

        adopt_init_args(self, locals())


class WriteHierarchicalModelSummaryTask(object):


    debug = False

    level_atoms_pdb = 'level_{:04d}-partitions.pdb'
    level_atoms_eff = 'level-definitions-level_{:04d}.eff'
    level_atoms_pymol_png_prefix = 'level_{:04d}-partitions'
    level_partitions_png = 'all-level-partitions-chain_{}.png'

    pymol_script_py = 'pymol_script.py'

    show_file_dict = show_file_dict

    def __init__(self,
        output_directory,
        master_phil,
        pymol_images = None,
        ):
        adopt_init_args(self, locals())

    def filepath(self, filename):
        return os.path.join(self.output_directory, filename)

    def run(self,
        reference_hierarchy,
        level_group_array,
        level_group_selection_strings,
        level_labels,
        overall_atom_mask,
        plotting_object,
        ):
        """Write out the composition of the hierarchical model"""

        logger.subheading('Writing summary of hierarchical model')

        # Object to generate structures with different b-factors etc.
        from pandemic.adp.hierarchy.utils import StructureFactory
        s_fac = StructureFactory(master_h=reference_hierarchy)

        output_files = collections.OrderedDict()

        of, level_hierarchies = self.write_level_atoms_structures(
            level_group_array = level_group_array,
            level_labels = level_labels,
            overall_atom_mask = overall_atom_mask,
            structure_factory = s_fac,
            )
        self.show_file_dict(of)
        output_files.update(of)

        of = self.make_level_atoms_eff_files(
            level_labels = level_labels,
            level_group_selection_strings = level_group_selection_strings,
            )
        self.show_file_dict(of)
        output_files.update(of)

        of = self.make_level_atoms_plots(
            level_hierarchies = level_hierarchies,
            level_labels = level_labels,
            structure_factory = s_fac,
            plotting_object = plotting_object,
            )
        self.show_file_dict(of)
        output_files.update(of)

        # Generate images of each chain of each level coloured by group
        if self.pymol_images is not None:
            of = self.make_pymol_images_of_level_atoms(
                level_atoms_pdb = output_files.get('level_atoms_pdb'),
                level_group_selection_strings = level_group_selection_strings,
                level_labels = level_labels,
                structure_factory = s_fac,
                )
            self.show_file_dict(of)
            output_files.update(of)

        of = {'pymol_script' : self.make_pymol_script(file_dict=output_files)}
        self.show_file_dict(of)
        output_files.update(of)

        self.result = HierarchySummaryOutput(
            output_files = output_files,
            )

        return self.result

    def write_level_atoms_structures(self,
        level_group_array,
        level_labels,
        overall_atom_mask,
        structure_factory,
        ):

        file_dict = collections.OrderedDict()

        # Generate hierarchy for each level with groups as b-factors
        level_hierarchies = []
        for i_level, g_vals in enumerate(level_group_array):
            level_lab = level_labels[i_level]
            all_values = -1 * numpy.ones_like(overall_atom_mask)
            all_values[overall_atom_mask] = g_vals # initially indexed from zero!
            m_h = structure_factory.custom_copy(uij=None, iso=all_values, mask=None, blank_copy=True)
            m_f = self.filepath(self.level_atoms_pdb.format(i_level+1))
            m_h.write_pdb_file(m_f)
            file_dict[level_lab] = m_f
            # Append to list for plotting
            level_hierarchies.append(m_h)

        return {'level_atoms_pdb' : file_dict}, level_hierarchies

    def make_level_atoms_eff_files(self,
        level_labels,
        level_group_selection_strings,
        ):

        from pandemic.adp.hierarchy.custom_levels import MakeNewCustomHierarchyEffFilesFromIndices
        write_levels_function = MakeNewCustomHierarchyEffFilesFromIndices(
            selection_strings_dict = dict(zip(level_labels, level_group_selection_strings)),
            master_phil = self.master_phil,
            custom_level_scope_name = 'model.custom_level',
            )

        file_dict = collections.OrderedDict()

        for i_l, (l, ls) in enumerate(zip(level_labels, level_group_selection_strings)):

            filename = os.path.join(self.output_directory, self.level_atoms_eff.format(i_l+1))

            write_levels_function(
                indices_hierarchy = [[[i] for i in range(len(ls))]],
                input_level_names = [l],
                output_level_names = [l],
                output_filename = filename,
                )

            file_dict[l] = filename

        return {'level_atoms_eff' : file_dict}

    def make_level_atoms_plots(self,
        level_hierarchies,
        level_labels,
        structure_factory,
        plotting_object,
        ):

        file_dict = collections.OrderedDict()

        # Write hierarchy plot for each chain
        b_h = structure_factory.blank_copy()
        b_c = b_h.atom_selection_cache()
        for c in b_h.chains():
            chain_sel = b_c.selection('chain {}'.format(c.id))
            hierarchies = [h.select(chain_sel, copy_atoms=True) for h in level_hierarchies]
            # Skip if no partitions in this chain
            #if (numpy.array([h.atoms().extract_b() for h in hierarchies]) == -1).all():
            #    continue
            filename = self.filepath(self.level_partitions_png.format(c.id))
            plotting_object.level_plots(
                filename=filename,
                hierarchies=hierarchies,
                labels = ['Level {}\n({})'.format(i_l+1, l) for i_l, l in enumerate(level_labels)],
                title='chain {}'.format(c.id),
                )
            file_dict[c.id] = filename

        return {'level_partitions_png' : file_dict}

    def make_pymol_images_of_level_atoms(self,
        level_atoms_pdb,
        level_group_selection_strings,
        level_labels,
        structure_factory,
        ):

        # Output file dict
        file_dict = collections.OrderedDict()

        # Record any failed/missing files
        missing_files = []

        from giant.structure.formatting import PymolSelection
        for level_lab, structure_filename in level_atoms_pdb.items():
            i_level = level_labels.index(level_lab)
            # Images for each chain (of the partitions) - coloured by group
            f_prefix = self.filepath(self.level_atoms_pymol_png_prefix.format(i_level+1))
            # Choose the style based on whether interested in atoms or larger regions
            styles = 'cartoon' if level_lab in ['chain','groups','sec. struct.'] else 'lines+spheres'
            # Create selections for each group in each level
            pymol_selections = translate_phenix_selections_to_pymol_selections_simple(level_group_selection_strings[i_level])
            # If returned selection is none, create atom-by-atom pymol selection
            blank_h = structure_factory.blank_copy()
            cache_h = blank_h.atom_selection_cache()
            selections = [s1 if (s1 is not None) else PymolSelection.join_or([PymolSelection.format(a) for a in blank_h.atoms().select(cache_h.selection(s2))]) for s1,s2 in zip(pymol_selections, level_group_selection_strings[i_level])]
            # Create image of different selections
            of = auto_chain_images(
                structure_filename = structure_filename,
                output_prefix = f_prefix,
                style = styles,
                het_style = 'lines+spheres',
                colours = ['red','green','blue'],
                colour_selections = selections,
                settings = [
                    ('cartoon_oval_length', 0.5),
                    ('cartoon_discrete_colors', 'on'),
                    ('sphere_scale', 0.25),
                    ],
                width=1000, height=750,
            )
            # Check output
            if (not of):
                logger.warning('no images have been generated: {}...'.format(f_prefix))
            else:
                # Append missing files
                [missing_files.append(v) for v in of.values() if not os.path.exists(v)]
                # Store in output dictionary
                file_dict[level_lab] = of

        if missing_files:
            logger.warning('\n'.join(['image does not exist: {}'.format(v) for v in missing_files]))

        return {'level_atoms_pymol_png' : file_dict}

    def make_pymol_script(self, file_dict):

        from giant.pymol_utils import PythonScript

        s = PythonScript(pretty_models=False, pretty_maps=False)

        s.change_into_directory_maybe(path=os.path.abspath(self.output_directory))

        for f in file_dict.get('level_atoms_pdb',{}).values():
            if f is None: continue
            obj = os.path.basename(f)
            s.load_pdb(
                f_name = os.path.relpath(os.path.abspath(f), start=self.output_directory),
                obj = os.path.basename(f),
                )
            s.colour(obj=obj, colour='grey')
            s.custom('spectrum', expression='b%10', palette="blue_white_green", selection="{} and (b>-1)".format(obj))

        s.show_as(obj='all', style='spheres')
        s.show(obj='all', style='sticks')
        s.set('sphere_scale', 0.25)
        s.set('grid_mode', 1)
        s.orient(obj='all')

        filename = self.filepath(self.pymol_script_py)
        s.write_script(filename)

        return filename


