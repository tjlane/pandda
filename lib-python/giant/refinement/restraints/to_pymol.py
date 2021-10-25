import giant.logs as lg
logger = lg.getLogger(__name__)

import copy

from giant.pymol_utils import (
    PythonScript,
    shapes,
    )

from giant.structure.formatting import (
    PhenixSelection,
    )


class WriteRestraintsPymolScript(object):

    def __init__(self, filepath=None, base_python_script=None):

        self.filepath = filepath
        self.base_python_script = base_python_script

    def __call__(self, restraints_collection, hierarchy, filepath=None):

        logger.subheading('Writing restraints as pymol script')

        if filepath is None:
            filepath = self.filepath

        python_script = self.make(
            restraints_collection = restraints_collection,
            hierarchy = hierarchy,
            )

        if filepath is not None:

            python_script.write_script(
                f_name = str(filepath),
                )

            logger('Pymol script written to {}'.format(str(filepath)))

        return python_script

    def initialise(self):

        python_script = PythonScript(
            pretty_models = False,
            pretty_maps = False,
            )

        if self.base_python_script is not None:
            python_script.lines.extend(
                copy.deepcopy(self.base_python_script.lines)
                )

        return python_script

    def make(self, restraints_collection, hierarchy):

        python_script = self.initialise()

        self.make_distance_restraint_objects(
            python_script = python_script,
            restraints_collection = restraints_collection,
            hierarchy = hierarchy,
            )

        self.make_occupancy_restraint_objects(
            python_script = python_script,
            restraints_collection = restraints_collection,
            hierarchy = hierarchy,
            )

        return python_script

    def make_distance_restraint_objects(self, python_script, restraints_collection, hierarchy):

        h_cache = hierarchy.atom_selection_cache()

        r_shapes = {}

        for r in restraints_collection.distance_restraints:

            if r.length <= 1e-3:
                continue

            color = [0.,1.,1.]

            atom1 = hierarchy.select(h_cache.selection(PhenixSelection.format(r.atom1))).only_atom()
            atom2 = hierarchy.select(h_cache.selection(PhenixSelection.format(r.atom2))).only_atom()

            if r.atom1['altloc'] != '':
                r_hash = 'alt_'+str(r.atom1['altloc'])
            elif r.atom2['altloc'] != '':
                r_hash = 'alt_'+str(r.atom2['altloc'])
            else:
                r_hash = 'main'

            r_shapes.setdefault(r_hash,[]).extend([
                shapes.Sphere(
                    centre = tuple(atom1.xyz),
                    radius = 0.2,
                    color = color,
                    ),
                shapes.Cylinder(
                    start = tuple(atom1.xyz),
                    end = tuple(atom2.xyz),
                    radius = 0.1,
                    colors = [color, color],
                    ),
                shapes.Sphere(
                    centre = tuple(atom2.xyz),
                    radius = 0.2,
                    color = color,
                    ),
                ])

        for k, k_shapes in sorted(r_shapes.items()):
            python_script.add_shapes(
                k_shapes,
                obj = 'distance_restraints_'+str(k),
                )

        return None

    def make_occupancy_restraint_objects(self, python_script, restraints_collection, hierarchy):

        h_cache = hierarchy.atom_selection_cache()
        h_atoms = hierarchy.atoms()

        for i_r, r in enumerate(restraints_collection.occupancy_restraints):

            r_label = "OccupancyRestraint{i:03d}".format(i=i_r+1)

            for i_g, g in enumerate(r.occupancy_groups):

                g_label = "{l}_Group{i:02d}".format(l=r_label, i=i_g+1)

                selection_string = PhenixSelection.join_or(
                    [PhenixSelection.format(o) for o in g.objects]
                    )

                atom_selection = h_cache.selection(selection_string)

                g_shapes = [
                    shapes.Sphere(
                        centre = tuple(a.xyz),
                        radius = 0.1,
                        color = [1.0,1.0,1.0],
                        )
                    for a in h_atoms.select(atom_selection)
                ]

                python_script.add_shapes(
                    g_shapes,
                    obj = g_label,
                    )

        return None








