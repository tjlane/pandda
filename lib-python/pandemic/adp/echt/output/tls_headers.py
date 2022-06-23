from libtbx import adopt_init_args
import operator
from functools import reduce


class MultiModelTLSHeaderFactory(object):


    def __init__(self,
        tls_group_array,
        tls_selection_strings,
        tls_objects,
        overall_selection,
        ):

        adopt_init_args(self, locals())

        # Create selection strings for everything not in a tls group
        null_selection_strings = []
        for level_strings in tls_selection_strings:
            null_selection = "not ({})".format(
                " or ".join(["({})".format(s) for s in level_strings])
            )
            null_selection_strings.append(null_selection)
        self.null_selection_strings = null_selection_strings

    def __call__(self,
        i_levels = None,
        ):

        import mmtbx.tls.tools
        if i_levels is None:
            i_levels = list(range(len(self.tls_objects)))

        assert isinstance(i_levels, list)
        assert set(map(type, i_levels)) == {int}
        i_levels = sorted(i_levels)

        # Unique combinations of the tls groups across levels
        unique_tls_combinations = sorted(set(map(tuple,self.tls_group_array[i_levels].T.tolist())))

        # List of tlso objects for each dataset
        group_tlsos = []
        group_strings = []

        # Iterate through the different UNIQUE grouping combinations
        for i_groups in unique_tls_combinations:

            # Ignore null groups
            if i_groups.count(-1) == len(i_groups):
                continue

            # Level-group pairs
            i_level_groups = list(zip(i_levels, i_groups))

            # Combine the selection strings for this combination
            group_selection_strings = []
            # Overall selection string
            if (self.overall_selection is not None):
                group_selection_strings.append(self.overall_selection)
            # For each group, or the null string
            for i_l, i_g in i_level_groups:
                if (i_g != -1):
                    group_selection_strings.append(self.tls_selection_strings[i_l][i_g])
                else:
                    group_selection_strings.append(self.null_selection_strings[i_l])
            combined_selection_string = '('+') and ('.join(group_selection_strings)+')'

            # Extract the tls parameters for each level
            tls_objects = [self.tls_objects[i_l][i_g] for i_l, i_g in i_level_groups if (i_g != -1)]
            # Expand out the tls objects for each level
            tls_parameters_exp = [[ma.expand() for ma in obj.tls_parameters] for obj in tls_objects]
            # Calculate the total for each dataset for each level
            tls_matrices_levels = [[reduce(operator.add, ml) for ml in zip(*mll)] for mll in tls_parameters_exp]
            # Sum matrices over levels
            tls_matrices = [reduce(operator.add, ms) for ms in zip(*tls_matrices_levels)]

            # Extract origins for these groups and check they're all the same (!)
            all_group_origins = [obj.origins for obj in tls_objects]
            for o in all_group_origins:
                assert o.as_double().all_eq(all_group_origins[0].as_double())
            tls_origins = all_group_origins[0]
            assert len(tls_matrices) == len(tls_origins)

            # Create tlso objects
            tlsos = [mmtbx.tls.tools.tlso(t=tuple(m.get('T')), l=tuple(m.get('L')), s=tuple(m.get('S')), origin=tuple(o)) for m,o in zip(tls_matrices, tls_origins)]

            group_tlsos.append(tlsos)
            group_strings.append(combined_selection_string)

        # Create header for each dataset
        from giant.logs import ListStream
        dataset_headers = []
        for tlsos in zip(*group_tlsos):
            # Store output in list stream
            l = ListStream()
            # Get the remark records
            mmtbx.tls.tools.remark_3_tls(tlsos=tlsos, selection_strings=group_strings, out=l)
            # And finally append
            dataset_headers.append(str(l))

        return dataset_headers


