import os
import numpy
from giant.pymol_utils import PythonScript, shapes


class base_pymol_script(object):


    def __init__(self,
        tls_groups,
        connectivity_matrix,
        model_structure,
        output_filename,
        ):

        self.centroids = {}

        ##################################################

        self.structure_obj = 'level_uijs'

        self.n_groups = len(tls_groups)

        self.group_centres = 'tls-group-centres'
        self.group_atoms = 'tls-group-atoms'

        self.max_radius = 0.5
        self.atom_size = self.max_radius / 5.

        self.bg_colour  = [0.6,0.6,1.0]

        self.output_filename = output_filename
        self.output_directory = os.path.dirname(os.path.abspath(output_filename))

        self.script = self._build(
            tls_groups = tls_groups,
            connectivity_matrix = connectivity_matrix,
            model_structure = model_structure,
            )

    def _build(self,
        tls_groups,
        connectivity_matrix,
        model_structure,
        ):

        ##################################################

        s = PythonScript(pretty_models=False, pretty_maps=False)

        # Create
        s.change_into_directory_maybe(path=self.output_directory)
        s.load_pdb(
            f_name = os.path.relpath(os.path.abspath(model_structure), start=self.output_directory),
            obj = self.structure_obj,
            )

        s.show_as(obj=self.structure_obj, style='lines')
        s.show(obj=self.structure_obj, style='ellipsoids')
        #s.colour(obj=self.structure_obj, colour="grey90")
        s.custom('spectrum', expression='b', selection=self.structure_obj)
        s.disable(obj=self.structure_obj)

        ##################################################

        at_shapes_sph = [shapes.Alpha(1.0)] # everything that follows has this applied
        at_shapes_cyl = [shapes.Alpha(0.5)]

        gp_shapes = [shapes.Alpha(1.0)]
        gp_atoms = []

        # Create shapes for each atom and join into groups
        for i in range(self.n_groups):
            ig = tls_groups[i]
            coords = ig.coordinates[0:ig.n_atoms] # select only first dataset
            centroid = coords.mean()

            # Contents of groups
            for xyz in coords:
                at_shapes_sph.append(shapes.Sphere(
                    centre=xyz,
                    radius=self.atom_size,
                    ))
                at_shapes_cyl.append(shapes.Cylinder(
                    start=centroid, end=xyz, colors=[self.bg_colour,self.bg_colour],
                    radius=self.atom_size/2.,
                    ))
            # Centres of groups
            gp_shapes.append(shapes.Sphere(
                centre=centroid,
                radius=1.25*self.max_radius,
                color=self.bg_colour,
                ))
            gp_atoms.append(shapes.Pseudoatom(
                pos=centroid,
                label='{}'.format(i+1),
                ))
            # Store for later
            self.centroids[i] = centroid

        # Add cgo objects
        at_shapes = at_shapes_sph + at_shapes_cyl
        s.add_shapes(at_shapes, obj=self.group_atoms)

        # Add group centres
        s.add_shapes(gp_shapes, obj=self.group_centres)
        # Add labels to group centres
        for o in gp_atoms:
            s.add_generic_object(o, obj=self.group_centres+'-labels')
        #s.show_as(style='spheres', obj=self.group_centres)
        #s.colour(colour=self.bg_colour, obj=self.group_centres)
        s.set('label_position', (0.0, 4.0*self.max_radius, 0.0))
        s.set('label_size', 28) # negative -> angstroms
        s.label(selection=self.group_centres, expression="label")

        ##################################################

        # Connections between neighbouring groups
        con_shapes = []
        for i in range(self.n_groups):
            for j in range(self.n_groups):
                # Only need to go one-way
                if i == j:
                    break
                # Only process if proximate
                if not connectivity_matrix[i,j]:
                    continue
                # Create links to show similarity
                cyl = shapes.Cylinder(
                    start=self.centroids[i], end=self.centroids[j], colors=[self.bg_colour,self.bg_colour],
                    radius=0.5*self.max_radius,
                    )
                con_shapes.append(cyl)

        s.add_shapes(con_shapes, obj='tls-group-neighbours')
        s.disable(obj='tls-group-neighbours')

        ##################################################

        return s

    def add_comparison_matrix(self,
        comparison_matrix,
        colour_function,
        radius_function,
        obj_name,
        min_max_values = None,
        ):

        s = self.script

        # Create links to show similarity
        for i in range(self.n_groups):
            com_shapes = []
            for j in range(self.n_groups):
                # Extract tls similarity
                v = comparison_matrix[i,j]
                # Map similarity to colour value (e.g. 0-1 -> red-green)
                col = colour_function(v)
                rad = radius_function(v)
                # Create links to show similarity
                cyl = shapes.Cylinder(
                    start=self.centroids[i],
                    end=self.centroids[j],
                    colors=[col,col],
                    radius=rad,
                    )
                com_shapes.append(cyl)
            s.add_shapes(com_shapes, obj=obj_name, state=i+1)

        # Find min/max for colour bar
        if min_max_values is None:
            min_v = 0.0
            max_v = comparison_matrix.max()
        else:
            min_v, max_v = min_max_values

        # Spacing of colour samples
        del_v = (max_v - min_v) / 10.0

        if del_v > 0.0:
            color_range = numpy.arange(min_v, max_v+del_v, del_v)
        else:
            color_range = [min_v, max_v]

        s.ramp_new(
            obj=obj_name+'-scalebar',
            mol_or_map_obj=self.structure_obj,
            range = str([min_v, max_v]),
            color = str([list(colour_function(v)) for v in color_range]),
            )

        s.disable(obj=obj_name)
        s.disable(obj=obj_name+'-scalebar')

    def write(self):

        s = self.script

        s.rebuild()
        s.orient(obj='all')
        s.write_script(self.output_filename)


