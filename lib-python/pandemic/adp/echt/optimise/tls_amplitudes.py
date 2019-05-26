import copy
import numpy
from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure
from scitbx.array_family import flex


class ValidateTLSGroup_Amplitudes:


    def __init__(self,
            tls_amplitude_tol = -1,
            ):
        adopt_init_args(self, locals())

    def __call__(self,
            tls_amplitudes,
            ):
        return (tls_amplitudes < -self.tls_amplitude_tol).all_eq(False)


class OptimiseTLSGroup_Amplitudes_AdoptValues:


    def __init__(self,
            i_dataset,
            i_modes,
            ):
        n_parameters = len(i_modes)
        i_modes = flex.size_t(i_modes)
        i_dataset = flex.size_t([i_dataset])
        adopt_init_args(self, locals())

    def __call__(self,
            values,
            tls_parameters,
            ):
        assert len(values) == self.n_parameters
        for i_m, v in zip(self.i_modes, values):
            mdl = tls_parameters[i_m]
            mdl.amplitudes.set([v], self.i_dataset)


class OptimiseTLSGroup_Amplitudes_TargetEvaluator:


    start_target = 1e6

    def __init__(self,
            uij_target,
            target_function,
            other_target_functions,
            multi_dataset_tls_group,
            optimise_i_dataset,
            optimise_i_modes,
            tls_amplitude_tol,
            isotropic_mask = None,
            verbose = False,
            log = None,
            ):
        if log is None: log = Log()

        # Init values
        n_call = 0
        best_target = self.start_target
        best_target_terms = None

        # Function for inserting amplitudes into multi_dataset_tls_group
        adopt_values = OptimiseTLSGroup_Amplitudes_AdoptValues(
                i_dataset = optimise_i_dataset,
                i_modes = optimise_i_modes,
                )

        # Function to check e.g. amplitudes not negative
        validate = ValidateTLSGroup_Amplitudes(
                tls_amplitude_tol = tls_amplitude_tol,
                )

        adopt_init_args(self, locals(), exclude=['tls_amplitude_tol'])

        # Extract un-multiplied uijs from the input object
        d_sel_1d = flex.double(multi_dataset_tls_group.n_datasets, 0.0)
        d_sel_1d[optimise_i_dataset] = 1.0
        d_sel_2d = (d_sel_1d.matrix_outer_product(flex.double(multi_dataset_tls_group.n_atoms,1.0)) > 0.0).as_1d()
        base_uijs = [multi_dataset_tls_group.tls_parameters[i_m].matrices.uijs(
                sites_cart = flex.vec3_double(multi_dataset_tls_group.coordinates.select(d_sel_2d)),
                origin = multi_dataset_tls_group.origins[optimise_i_dataset],
                ) for i_m in optimise_i_modes]

        # Make some uijs isotropic if required
        if isotropic_mask is not None:
            base_uijs = [isotropic_mask(u) for u in base_uijs]

        # Convert to numpy arrays for convenience
        self.base_uijs = [u.as_1d().as_double().as_numpy_array() for u in base_uijs]

    def target(self,
            parameters,
            ):
        self.n_call += 1
        if self.verbose and (self.n_call%20==1): self.print_header_line(parameters)
        # Validate the new combination of parameters and return poor score if invalid
        if not self.validate(tls_amplitudes=parameters):
            if self.verbose: self.print_invalid_line(parameters)
            return 1.1 * self.best_target
        # Apply multipliers to each mode
        uij_modes = [(a*u).reshape((u.size//6,6)) for a, u in zip(parameters, self.base_uijs)]
        uij_fitted = flex.sym_mat3_double(numpy.sum(uij_modes, axis=0))
        uij_deltas = uij_fitted - self.uij_target
        # Calculate target function
        lsq = self.target_function(uij_deltas=uij_deltas)
        # Calculate full target and store if minimum
        target_terms = [lsq] + [f(lsq=lsq, uij_deltas=uij_deltas) for f in self.other_target_functions]
        target = sum(target_terms)
        if target < self.best_target:
            self.best_target = target
            self.best_target_terms = target_terms
        if self.verbose: self.print_current_line(parameters, target)
        return target

    def print_header_line(self, values):
        header = '[{}] -> ({:^12}, {:^10})'.format(', '.join(['{:>10}'.format('param') for r in values]), 'target', '(best)')
        line = '-'*len(header)
        self.log(line)
        self.log(header)
        self.log(line)

    def print_current_line(self, values, target):
        self.log('[{}] -> ({:12.9f}, {:10.7f})'.format(', '.join(['{:+10.7f}'.format(r) for r in values]), target, self.best_target))

    def print_invalid_line(self, values):
        self.log('[{}] -> ({:>12}, {:>10})'.format(', '.join(['{:+10.7f}'.format(r) for r in values]), 'UNPHYSICAL', ''))


class OptimiseTLSGroup_Amplitudes:


    evaluator_class = OptimiseTLSGroup_Amplitudes_TargetEvaluator

    def __init__(self,
            uij_target,
            uij_target_weights,
            simplex_delta,
            target_function,
            other_target_functions,
            convergence_tol,
            tls_amplitude_tol,
            isotropic_mask = None,
            verbose = False,
            log = None,
            ):

        # Convert to flex
        uij_target = [flex.sym_mat3_double(a) for a in uij_target]
        uij_target_weights = [flex.double(numpy.ascontiguousarray(a)) for a in uij_target_weights]

        if log is None: log = Log()
        adopt_init_args(self, locals(), exclude=['a'])

    def __call__(self,
            multi_dataset_tls_group,
            optimise_i_dataset,
            optimise_i_modes,
            ):

        # Extract target values for this dataset
        t_uij = self.uij_target[optimise_i_dataset]
        t_wgt = self.uij_target_weights[optimise_i_dataset]

        # Copy and set variables for the target functions
        tf = copy.deepcopy(self.target_function)
        otf = copy.deepcopy(self.other_target_functions)
        for f in [tf]+otf:
            f.set_uij_weights(t_wgt)

        # Initialise new instance for each call
        evaluator = self.evaluator_class(
                uij_target = t_uij,
                target_function = tf,
                other_target_functions = otf,
                multi_dataset_tls_group = multi_dataset_tls_group,
                optimise_i_dataset = optimise_i_dataset,
                optimise_i_modes = optimise_i_modes,
                tls_amplitude_tol = self.tls_amplitude_tol,
                isotropic_mask = self.isotropic_mask,
                verbose = self.verbose,
                log = self.log,
                )

        # Generate simplex
        starting_amplitudes = [multi_dataset_tls_group.tls_parameters[i_m].amplitudes[optimise_i_dataset] for i_m in optimise_i_modes]
        starting_simplex = [starting_amplitudes]
        for i in xrange(len(starting_amplitudes)):
            new_simplex_point = copy.copy(starting_amplitudes)
            new_simplex_point[i] = new_simplex_point[i] + self.simplex_delta
            starting_simplex.append(new_simplex_point)

        # Optimise these parameters
        from scitbx import simplex
        optimised = simplex.simplex_opt(dimension = len(starting_amplitudes),
                                        matrix    = map(flex.double, starting_simplex),
                                        evaluator = evaluator,
                                        tolerance = self.convergence_tol)

        # Adopt the solution!
        evaluator.adopt_values(
                values = optimised.get_solution(),
                tls_parameters = multi_dataset_tls_group.tls_parameters,
                )

        # Zero negative values to ensure validity
        multi_dataset_tls_group.tls_parameters.zero_negative_amplitudes()

        # Create results object
        results = group_args(
            n_iter = evaluator.n_call,
            target_function_value = evaluator.best_target,
            target_function_values = evaluator.best_target_terms,
            )

        return (multi_dataset_tls_group, results)


