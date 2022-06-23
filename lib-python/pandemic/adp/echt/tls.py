import numpy
from libtbx import adopt_init_args
from scitbx.array_family import flex


class EchtBFactorModel(object):


    def __init__(self,
            tls_objects,
            tls_selections,
            adp_values = None,
            tls_level_names = None,
            adp_level_name = 'atomic',
            tls_selection_strings = None,
            dataset_labels = None,
            ):
        n_tls_levels = len(tls_objects)
        n_levels = n_tls_levels + 1
        n_modes = tls_objects[0][0].n_modes
        n_datasets = tls_objects[0][0].n_datasets
        n_atoms = len(tls_selections[0][0])
        n_groups = 0 # will be populated later
        if adp_values is None:
            adp_values = flex.sym_mat3_double(n_atoms)
        else:
            adp_values = flex.sym_mat3_double(adp_values)
        if dataset_labels is None:
            dataset_labels = ['Dataset {}'.format(i+1) for i in range(n_datasets)]
        if tls_level_names is None:
            tls_level_names = ['Level {}'.format(i+1) for i in range(n_tls_levels)]
        all_level_names = tls_level_names + [adp_level_name]
        all_level_types = ['tls']*n_tls_levels + ['adp']
        adopt_init_args(self, locals())
        self.process_input()

    def process_input(self):
        assert self.n_tls_levels > 0
        assert self.n_datasets > 0
        assert self.n_atoms > 0
        assert len(self.tls_objects) == self.n_tls_levels
        assert len(self.tls_selections) == self.n_tls_levels
        assert len(self.tls_level_names) == self.n_tls_levels
        assert len(self.adp_values) == self.n_atoms
        assert len(self.all_level_names) == self.n_levels
        assert len(set(self.all_level_names)) == self.n_levels, 'Some level names are duplicated -- please rename level names so they are unique: {}'.format(self.all_level_names)
        assert len(self.dataset_labels) == self.n_datasets
        for objs, sels in zip(self.tls_objects, self.tls_selections):
            # Counter of number of groups in hierarchy
            self.n_groups += len(objs)
            assert len(objs) == len(sels)
            for o, s in zip(objs, sels):
                assert o.n_datasets == self.n_datasets
                assert o.n_modes == self.n_modes
                assert s.size == self.n_atoms
                assert s.sum() == o.n_atoms

    def uijs(self):
        """
        Returns the uij values for each level for each dataset
        Size of output array: (n_levels, n_datasets, n_atoms, 6)
        """
        uij_values = numpy.zeros(shape=(self.n_levels, self.n_datasets, self.n_atoms, 6))
        # Iterate through levels and extract tls values from each group in each level
        for i_level in range(self.n_tls_levels):
            # Iterate through groups in level
            for g_sel, g_obj in zip(self.tls_selections[i_level], self.tls_objects[i_level]):
                g_uijs = numpy.array(g_obj.uijs()).reshape(
                    (g_obj.n_datasets, g_obj.n_atoms, 6),
                    )
                uij_values[i_level][:, g_sel] = g_uijs
        # Store the fixed adps as the last level
        for i_dataset in range(self.n_datasets):
            uij_values[self.n_levels-1, i_dataset] = self.adp_values
        return uij_values

    def uijs_tls(self):
        """
        Returns the uij values for each tls level for each dataset.
        Size of output array: (n_tls_levels, n_modes, n_datasets, n_atoms, 6)
        """
        uij_values = numpy.zeros(shape=(self.n_tls_levels, self.n_modes, self.n_datasets, self.n_atoms, 6))
        # Iterate through levels
        for i_level in range(self.n_tls_levels):
            # Iterate through groups in level
            for g_sel, g_obj in zip(self.tls_selections[i_level], self.tls_objects[i_level]):
                g_uijs = numpy.array(g_obj.uijs_by_mode()).reshape(
                    (g_obj.n_modes, g_obj.n_datasets, g_obj.n_atoms, 6),
                    )
                uij_values[i_level][:, :, g_sel] = g_uijs
        return uij_values


class MultiDatasetTLSGroup(object):


    def __init__(self,
            index,
            label,
            tls_parameters,
            coordinates,
            origins,
            ):
        adopt_init_args(self, locals())
        self.convert_input_from_numpy_to_flex()
        self.n_modes = self.tls_parameters.size()
        self.n_datasets, self.n_atoms = self.coordinates.all()
        self.validate_input()

    def convert_input_from_numpy_to_flex(self):
        if hasattr(self.origins, 'shape'):
            self.origins = flex.vec3_double(self.origins)
        else:
            assert hasattr(self.origins, 'all'), 'input origins must be flex or numpy array'
        if hasattr(self.coordinates, 'shape'):
            shape = self.coordinates.shape
            self.coordinates = flex.vec3_double(self.coordinates.reshape(numpy.product(shape[:-1]),3))
            self.coordinates.reshape(flex.grid(shape[:-1]))
        else:
            assert hasattr(self.coordinates, 'all'), 'input coordinates must be flex or numpy array'

    def validate_input(self):
        assert self.n_modes > 0
        assert self.n_datasets > 0
        assert self.n_atoms > 0
        assert self.origins.nd() == 1
        assert self.coordinates.nd() == 2
        assert self.origins.size() == self.n_datasets
        for mode in self.tls_parameters:
            assert mode.amplitudes.size() == self.n_datasets

    def uijs(self):
        return self.tls_parameters.uijs(
                sites_carts = self.coordinates,
                origins = self.origins,
                )

    def uijs_by_mode(self):
        return [m.uijs(
                sites_carts = self.coordinates,
                origins = self.origins,
                ) for m in self.tls_parameters]

    def uijs_unmultiplied(self):
        """Extract amplitudes and uijs as separate objects (would normally by pairwise multiplied)"""
        return [(m.amplitudes.get(), [m.matrices.uijs(
                sites_cart = self.coordinates[i*self.n_atoms:(i+1)*self.n_atoms],
                origin = self.origins[i],
                ) for i in range(self.n_datasets)],
                ) for m in self.tls_parameters]

    def summary(self):
        """Print the number of parameters/input data"""
        bar = '=========>'
        s = bar+'\nTLS Group Fit Summary: {}\n'.format(self.label)+bar
        for i_tls in range(self.n_models):
            s += '\n> TLS model {}'.format(i_tls+1)
            mode = self.tls_parameters.get(index=i_tls)
            s += '\n\t' + mode.matrices.summary().replace('\n','\n\t')
            s += '\n\t' + mode.amplitudes.summary().replace('\n','\n\t')
        return s
