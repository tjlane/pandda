import tqdm
from libtbx import adopt_init_args
from libtbx.utils import Sorry, Failure
from scitbx.array_family import flex
from bamboo.common.logs import Log
from pandemic.adp.echt.optimise.targets import TargetTerm_UijLeastSquares, TargetTerm_UijTargetBarrier


class OptimiseTLSGroup:


    target_function_class = TargetTerm_UijLeastSquares
    barrier_function_class = TargetTerm_UijTargetBarrier

    std_component_optimisation_list = ["TLS"]
    null_component_optimisation_list = ["T","L","S","TL","LS","TS","TLS"]

    def __init__(self,
            convergence_tol,
            tolerances,
            eps_values,
            simplex_params,
            barrier_params=None,
            #n_cycles=5,
            n_cpus=1,
            verbose = False,
            log = None,
            ):

        if eps_values is None:
            eps_values = group_args(
                    tls_matrices_eps = -1,
                    tls_amplitudes_eps = -1,
                    )

        # Use selected target function
        target_function = self.target_function_class()

        # Construct additional target terms
        other_target_functions = []
        if barrier_params is not None:
            other_target_functions.append(
                    self.barrier_function_class(
                        barrier_height = barrier_params.barrier_height,
                        barrier_width = barrier_params.barrier_width,
                        barrier_offset = barrier_params.barrier_offset,
                        )
                    )

        #ls = LogStream(path=os.path.splitext(self.log.log_file().path)[0]+'-level{:04d}-group{:06d}.log'.format(self.index, i))
        #log = Log(stdout=True).add_output(ls).toggle(status=(self._n_groups==1))
        adopt_init_args(self, locals())

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

        from pandemic.adp.echt.optimise.tls_matrices import OptimiseTLSGroup_Matrices, TLSSimplexGenerator
        from pandemic.adp.echt.optimise.tls_amplitudes import OptimiseTLSGroup_Amplitudes

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
            convergence_tol = self.convergence_tol,
            tls_matrix_tol = self.tolerances.tls_matrix_tolerance,
            tls_matrix_eps = self.eps_values.tls_matrices_eps,
            isotropic_mask = uij_isotropic_mask,
            verbose = self.verbose,
            log = self.log,
            )

        # Create optimisation task for amplitude optimisation
        # optimise_amplitudes = OptimiseTLSGroup_Amplitudes(
        #     uij_target = uij_target,
        #     uij_target_weights = uij_target_weights,
        #     simplex_delta = self.simplex_params.amplitude_delta,
        #     target_function = self.target_function,
        #     other_target_functions = self.other_target_functions,
        #     convergence_tol = self.convergence_tol,
        #     tls_amplitude_tol = self.tolerances.tls_amplitude_tolerance,
        #     isotropic_mask = uij_isotropic_mask,
        #     verbose = self.verbose,
        #     log = self.log,
        #     )

        # Extract loop parameters
        i_modes = range(multi_dataset_tls_group.tls_parameters.size())

        # Do cycles of alternating optimisation
        for i_mode in i_modes:

            # Extract mode object
            mode = multi_dataset_tls_group.tls_parameters[i_mode]

            # Skip if null
            if mode.is_null(
                    matrices_tolerance = self.eps_values.tls_matrices_eps,
                    amplitudes_tolerance = self.eps_values.tls_amplitudes_eps,
                    ):
                self.reset_mode(mode)
                continue

            components_list = self.null_component_optimisation_list

            # Iterate through component set
            for optimise_components in components_list:
                # Run matrix optimisation
                multi_dataset_tls_group, _ = optimise_matrices(
                        multi_dataset_tls_group = multi_dataset_tls_group,
                        optimise_components = optimise_components,
                        optimise_i_mode = i_mode,
                        )

            # If has refined to zero-values, reset.
            if mode.is_null(
                    matrices_tolerance = self.eps_values.tls_matrices_eps,
                    amplitudes_tolerance = self.eps_values.tls_amplitudes_eps,
                    ):
                self.reset_mode(mode)
                continue

            # Check that mode is valid after normalisation
            if not mode.matrices.is_valid():

                # Attempt to fix by reoptimising
                multi_dataset_tls_group, _ = optimise_matrices(
                        multi_dataset_tls_group = multi_dataset_tls_group,
                        optimise_components = "TLS",
                        optimise_i_mode = i_mode,
                        )

                # If still not valid, reset and re-optimise
                if not mode.matrices.is_valid():
                    mode.matrices.reset()
                    multi_dataset_tls_group, _ = optimise_matrices(
                            multi_dataset_tls_group = multi_dataset_tls_group,
                            optimise_components = "TLS",
                            optimise_i_mode = i_mode,
                            )

                # Reset mode
                if not mode.matrices.is_valid():
                    self.reset_mode(mode)
                    continue

            # Normalise the resulting matrices and amplitudes
            mode.normalise_by_matrices(
                    sites_carts=multi_dataset_tls_group.coordinates,
                    origins=multi_dataset_tls_group.origins,
                    target=1.0,
                    )

            # Check that mode is valid after normalisation
            if not mode.matrices.is_valid():

                # Attempt to fix by reoptimising
                multi_dataset_tls_group, _ = optimise_matrices(
                        multi_dataset_tls_group = multi_dataset_tls_group,
                        optimise_components = "TLS",
                        optimise_i_mode = i_mode,
                        )

                # If still not valid, reset and re-optimise
                if not mode.matrices.is_valid():
                    mode.matrices.reset()
                    multi_dataset_tls_group, _ = optimise_matrices(
                            multi_dataset_tls_group = multi_dataset_tls_group,
                            optimise_components = "TLS",
                            optimise_i_mode = i_mode,
                            )

                # Reset mode
                if not mode.matrices.is_valid():
                    self.reset_mode(mode)
                    continue
                    # raise Failure('Failed to fix core model. \n{}'.format(mode.summary()))

        # If everything is null, stop optimisation
        if multi_dataset_tls_group.tls_parameters.is_null(
                matrices_tolerance = self.eps_values.tls_matrices_eps,
                amplitudes_tolerance = self.eps_values.tls_amplitudes_eps,
                ):
            for i_mode, mode in enumerate(multi_dataset_tls_group.tls_parameters):
                self.reset_mode(mode)
                if i_mode != 0:
                    mode.matrices.set(values=(0.)*21, component_string='TLS')

        # TODO
        optimisation_info = None

        # Only return the minimum to reduce pickle overhead
        return (multi_dataset_tls_group, optimisation_info)

    @staticmethod
    def reset_mode(mode):
        mode.matrices.reset()
        mode.matrices.set(values=(1.,1.,1.,0.,0.,0.), component_string='T')
        mode.amplitudes.zero_values()


class OptimiseTLSLevel:


    def __init__(self,
            optimisation_function,
            n_cpus = 1,
            log = None,
            ):
        if log is None: log = Log()
        from pandemic.adp.parallel import MultiProcessWrapper, RunParallelWithProgressBarUnordered
        optimisation_wrapper = MultiProcessWrapper(function=optimisation_function)
        run_parallel = RunParallelWithProgressBarUnordered(n_cpus)
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
            ) for i_grp, (grp, sel) in enumerate(zip(tls_groups, tls_selections))]

        if (self.n_cpus == 1):
            #results = [self.optimisation_function(**a) for a in args]
            results = []
            pbar = tqdm.tqdm(total=len(args), ncols=100)
            for a in args: 
                results.append(self.optimisation_function(**a))
                pbar.update(1)
            pbar.close()
        else:
            # update chunksize for parallelisation
            self.run_parallel.chunksize = min(100, max(1, int(len(args)/float(self.run_parallel.n_cpus))))

            [a.update({'sort_value':i}) for i,a in enumerate(args)]
            results = self.run_parallel(function=self.optimisation_wrapper, arg_list=args)
            results = [r[1] for r in sorted(results, key=lambda x: x[0])]

        errors = []
        for r in results:
            if isinstance(r, str):
                self.log.bar()
                self.log(r)
                errors.append(r)

        if errors:
            self.log.bar()
            raise Failure('{} errors raised during optimisation (above)'.format(len(errors)))

        new_tls_groups = [r[0] for r in results]
        output_group_labels = [g.label for g in new_tls_groups]
        assert input_group_labels == output_group_labels

        return new_tls_groups



