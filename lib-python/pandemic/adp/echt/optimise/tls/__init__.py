import giant.logs as lg
logger = lg.getLogger(__name__)

import tqdm
from libtbx import adopt_init_args
from libtbx.utils import Sorry, Failure
from scitbx.array_family import flex
from pandemic.adp.echt.optimise.targets import TargetTerm_UijLeastSquares


class OptimiseTLSGroup(object):


    target_function_class = TargetTerm_UijLeastSquares

    component_optimisation_list = ["T","L","S","TL","LS","TS","TLS"]

    def __init__(self,
        simplex_params,
        convergence_tolerance,
        eps_values,
        ):

        # if eps_values is None:
        #     eps_values = group_args(
        #         matrix_eps = -1,
        #         amplitude_eps = -1,
        #         )

        # Use selected target function
        target_function = self.target_function_class()

        # Construct additional target terms
        other_target_functions = []

        adopt_init_args(self, locals())

        self.validate()

    def validate(self):

        for f in ['vibration_delta', 'libration_delta']:
            assert hasattr(self.simplex_params, f)

        for f in ['matrix_eps', 'amplitude_eps']:
            assert hasattr(self.eps_values, f)

    def __call__(self,
        multi_dataset_tls_group,
        uij_target,
        uij_target_weights,
        uij_isotropic_mask,
        uij_optimisation_mask,
        ):

        n_datasets, n_atoms = uij_target_weights.shape

        assert uij_target.shape == (n_datasets, n_atoms, 6)
        assert uij_target_weights.shape == (n_datasets, n_atoms)

        from pandemic.adp.echt.optimise.tls.matrices import OptimiseTLSGroup_Matrices, TLSSimplexGenerator

        # Create optimisation task for matrix optimisation
        optimise_matrices = OptimiseTLSGroup_Matrices(
            uij_target = uij_target,
            uij_target_weights = uij_target_weights,
            simplex_generator = TLSSimplexGenerator(
                vibration_delta = self.simplex_params.vibration_delta,
                libration_delta = self.simplex_params.libration_delta,
                ),
            target_function = self.target_function,
            other_target_functions = self.other_target_functions,
            convergence_tolerance = self.convergence_tolerance,
            uij_isotropic_mask = uij_isotropic_mask,
            )

        # Do cycles of alternating optimisation
        for i_mode in range(multi_dataset_tls_group.tls_parameters.size()):

            # Extract mode object
            mode = multi_dataset_tls_group.tls_parameters[i_mode]

            # Extract the input amplitudes to apply at the end
            initial_amplitudes = mode.amplitudes.get().deep_copy()

            # Skip if amplitudes are exactly zero
            if mode.amplitudes.get().all_eq(0.0):
                continue

            # Reset at beginning if matrices are null or invalid
            if not mode.matrices.is_valid():
                # Just set the matrices to zero for following optimisation
                mode.matrices.reset()

            # Iterate through component set
            for cpts in self.component_optimisation_list:
                # Run matrix optimisation
                multi_dataset_tls_group, _ = optimise_matrices(
                    multi_dataset_tls_group = multi_dataset_tls_group,
                    optimise_components = cpts,
                    optimise_i_mode = i_mode,
                    )

            # # If has refined to zero-values, reset.
            # if mode.is_null(
            #         matrices_tolerance = self.eps_values.matrix_eps,
            #         amplitudes_tolerance = self.eps_values.amplitude_eps,
            #         ):
            #     self.default_matrix(mode)
            #     continue

            # Normalise the resulting matrices and amplitudes
            if mode.matrices.is_valid():
                ret = mode.normalise_by_matrices(
                    sites_carts = multi_dataset_tls_group.coordinates,
                    origins = multi_dataset_tls_group.origins,
                    target = 1.0,
                    )
                # matrices produce negligible uijs -- continue
                if ret == -1:
                    self.zero_amplitudes(mode)
                    self.default_matrix(mode)
                    continue

            # Check that mode is valid after normalisation
            if not mode.matrices.is_valid():

                # Attempt to fix by reoptimising
                for cpts in self.component_optimisation_list:
                    multi_dataset_tls_group, _ = optimise_matrices(
                        multi_dataset_tls_group = multi_dataset_tls_group,
                        optimise_components = cpts,
                        optimise_i_mode = i_mode,
                        )
                    # # Break as soon as it becomes valid
                    # if mode.matrices.is_valid():
                    #     break

                # If still not valid, reset and re-optimise (now with the amplitudes at the "correct" scale)
                if not mode.matrices.is_valid():
                    mode.matrices.reset()
                    for cpts in self.component_optimisation_list:
                        multi_dataset_tls_group, _ = optimise_matrices(
                            multi_dataset_tls_group = multi_dataset_tls_group,
                            optimise_components = cpts,
                            optimise_i_mode = i_mode,
                            )

            # Ensure normalised at the end
            if mode.matrices.is_valid():
                ret = mode.normalise_by_matrices(
                    sites_carts = multi_dataset_tls_group.coordinates,
                    origins = multi_dataset_tls_group.origins,
                    target = 1.0,
                    )
                # matrices produce negligible uijs -- continue
                if ret == -1:
                    self.zero_amplitudes(mode)
                    self.default_matrix(mode)

            # # # # # # # # # # # #
            # THIS IS "GIVING UP"
            #
            # At the end, if nothing else, reset everything, to ensure that something valid is returned...
            if not mode.matrices.is_valid():
                self.zero_amplitudes(mode)
                self.default_matrix(mode)
                continue
            #
            # THIS IS "GIVING UP"
            # # # # # # # # # # # #

            # At the end of the cycle, apply the initial amplitudes
            mode.amplitudes.set(initial_amplitudes)

        # If everything is null, reset everything (first mode non-zero T, otherwise all zero)
        # if multi_dataset_tls_group.tls_parameters.is_null(
        #         matrices_tolerance = self.eps_values.matrix_eps,
        #         amplitudes_tolerance = self.eps_values.amplitude_eps,
        #         ):
        #     for i_mode, mode in enumerate(multi_dataset_tls_group.tls_parameters):
        #         self.default_matrix(mode)
        #     # Give non-zero matrices to first mode
        #     self.default_matrix(multi_dataset_tls_group.tls_parameters[0])

        # TODO
        optimisation_info = None

        # Only return the minimum to reduce pickle overhead
        return multi_dataset_tls_group

    @staticmethod
    def default_matrix(mode):
        mode.matrices.reset()
        mode.matrices.set(values=(1.,1.,1.,0.,0.,0.), component_string='T')

    @staticmethod
    def zero_amplitudes(mode):
        mode.amplitudes.zero_values()


class OptimiseTLSLevel(object):


    def __init__(self,
            optimisation_function,
            n_cpus = 1,
            ):
        from pandemic.adp.parallel import RunParallelWithProgressBarUnordered
        run_parallel = RunParallelWithProgressBarUnordered(
            function = optimisation_function,
            n_cpus = n_cpus,
            max_chunksize = 10,
            keep_processes_open = True,
        )
        adopt_init_args(self, locals())

    def __call__(self,
            tls_groups,
            tls_selections,
            uij_target,
            uij_target_weights,
            uij_isotropic_mask,
            uij_optimisation_mask = None,
            ):

        input_group_labels = [g.label for g in tls_groups]

        # Arg list
        args = [dict(
            multi_dataset_tls_group = grp,
            uij_target = uij_target[:,sel],
            uij_target_weights = uij_target_weights[:,sel],
            uij_isotropic_mask = (uij_isotropic_mask[sel] if uij_isotropic_mask else None),
            uij_optimisation_mask = uij_optimisation_mask,
            ) for (grp, sel) in zip(tls_groups, tls_selections)]

        if (self.n_cpus == 1):
            results = []
            pbar = tqdm.tqdm(total=len(args), ncols=100)
            for a in args:
                results.append(self.optimisation_function(**a))
                pbar.update(1)
            pbar.close()
        else:
            # Run with parallel wrapper
            results = self.run_parallel(arg_dicts=args)

        errors = []
        for r in results:
            if isinstance(r, str):
                logger.bar()
                logger(r)
                errors.append(r)

        if errors:
            logger.bar()
            raise Failure('{} errors raised during optimisation (above)'.format(len(errors)))

        new_tls_groups = results
        output_group_labels = [g.label for g in new_tls_groups]
        assert input_group_labels == output_group_labels

        return new_tls_groups
