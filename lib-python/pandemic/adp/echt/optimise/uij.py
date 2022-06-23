import giant.logs as lg
logger = lg.getLogger(__name__)

import copy, tqdm
import numpy as np
from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure
from scitbx.array_family import flex
from pandemic.adp.echt.optimise.targets import TargetTerm_UijLeastSquares


class UijSimplexGenerator(object):


    def __init__(self,
            uij_delta,
            ):
        adopt_init_args(self, locals())

    def __call__(self,
            uij_values,
            ):

        import numpy

        start = uij_values
        n = 6
        assert n == len(start)

        simplex_deltas = numpy.zeros((n+1,n))
        for i in range(n):
            simplex_deltas[i+1,i] = self.uij_delta

        start_simplex = numpy.repeat([start], len(start)+1, axis=0)
        assert start_simplex.shape == simplex_deltas.shape
        start_simplex += simplex_deltas
        return start_simplex

    def summary(self):
        s = 'Simplex optimisation deltas'
        s += '\nUij Delta:  {}'.format(self.uij_delta)
        return s


from pandemic.adp.echt.validate.uij import ValidateUijValues
class OptimiseUijValue_TargetEvaluator(object):


    start_target = 1e6

    def __init__(self,
            uij_target,
            target_function,
            other_target_functions,
            uij_tol,
            ):

        # Init values
        n_call = 0
        best_target = self.start_target
        best_target_terms = None

        # Function to check e.g. amplitudes not negative
        validate = ValidateUijValues(uij_tol = uij_tol)

        adopt_init_args(self, locals(), exclude=['uij_tol'])

    def target(self,
            parameters,
            ):
        self.n_call += 1
        if not self.validate(u=parameters):
            return 1.1 * abs(self.best_target)
        # Apply multipliers to each mode
        uij_fitted = flex.sym_mat3_double(self.uij_target.size(), tuple(parameters))
        uij_deltas = uij_fitted - self.uij_target
        # Calculate target function
        lsq = self.target_function(uij_deltas=uij_deltas)
        # Calculate full target and store if minimum
        target_terms = [lsq] + [f(lsq=lsq, uij_deltas=uij_deltas) for f in self.other_target_functions]
        target = sum(target_terms)
        if target < self.best_target:
            self.best_target = target
            self.best_target_terms = target_terms
        return target


class OptimiseUijValue(object):


    target_function_class = TargetTerm_UijLeastSquares
    evaluator_class = OptimiseUijValue_TargetEvaluator

    def __init__(self,
            simplex_params,
            convergence_tolerance,
            uij_eps,
            uij_tolerance,
            ):

        # Create target terms
        target_function = self.target_function_class()
        other_target_functions = []

        # Generate simplexes for optimisation
        simplex_generator = UijSimplexGenerator(uij_delta = simplex_params.uij_delta)

        # Default return value for small values
        eps_sphere = (uij_eps, uij_eps, uij_eps, 0., 0., 0.)

        adopt_init_args(self, locals())

    def __call__(self,
            uij_value,
            uij_target,
            uij_target_weights,
            uij_optimisation_mask,
            ):

        # Process inputs
        uij_target = flex.sym_mat3_double(uij_target)
        uij_target_weights = flex.double(uij_target_weights)

        # Make copies of functions to avoid contagion
        target_function = copy.deepcopy(self.target_function)
        other_target_functions = copy.deepcopy(self.other_target_functions)
        for f in [target_function]+other_target_functions:
            f.set_uij_weights(uij_target_weights)

        # Create evaluator
        evaluator = self.evaluator_class(
                uij_target = uij_target,
                target_function = target_function,
                other_target_functions = other_target_functions,
                uij_tol = self.uij_tolerance,
                )

        # Check if uij is valid, else start from zero
        if not evaluator.validate(u=uij_value):
            uij_value = (0.,)*6

        # If small values, don't optimise and return spherical uij of size eps
        if (flex.abs(flex.double(uij_value)) < self.uij_eps).all_eq(True):
            return self.eps_sphere

        # Generate simplex
        starting_simplex = self.simplex_generator(uij_value)

        # Optimise these parameters
        from scitbx import simplex
        optimised = simplex.simplex_opt(dimension = 6, # len(starting_simplex[0]),
                                        matrix    = list(map(flex.double, starting_simplex)),
                                        evaluator = evaluator,
                                        tolerance = self.convergence_tolerance)

        # # Create results object
        # results = group_args(
        #     n_iter = evaluator.n_call,
        #     target_function_value = evaluator.best_target,
        #     target_function_values = evaluator.best_target_terms,
        #     )

        return optimised.get_solution()


class OptimiseUijLevel(object):


    def __init__(self,
            optimisation_function,
            n_cpus = 1,
            ):
        from pandemic.adp.parallel import RunParallelWithProgressBarUnordered
        run_parallel = RunParallelWithProgressBarUnordered(
            function = optimisation_function,
            n_cpus = n_cpus,
            max_chunksize = 100,
            keep_processes_open = True,
        )
        adopt_init_args(self, locals())

    def __call__(self,
            uij_values,
            uij_target,
            uij_target_weights,
            uij_isotropic_mask = None,
            uij_optimisation_mask = None,
            ):

        if uij_optimisation_mask is not None:
            uij_target = uij_target[uij_optimisation_mask]
            uij_target_weights = uij_target_weights[uij_optimisation_mask]

        # Arg list
        args = [dict(
            uij_value = u,
            uij_target = uij_target[:,i_u],
            uij_target_weights = uij_target_weights[:,i_u],
            uij_optimisation_mask = uij_optimisation_mask,
            ) for i_u, u in enumerate(uij_values)]

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

        # Repackage optimised values
        new_adp_values = flex.sym_mat3_double(list(map(tuple,results)))

        # Make necessary values isotropic
        if uij_isotropic_mask is not None:
            new_adp_values = uij_isotropic_mask(new_adp_values)

        # Scale output to input 
        new_adp_values = self.scale_atomwise_to_input(
            output_adps = new_adp_values,
            input_adps = uij_values,
            )

        return new_adp_values

    def scale_atomwise_to_input(self,
        output_adps, 
        input_adps,
        eps = 1e-3,
        ):

        output_adps = np.array(output_adps)
        input_adps = np.array(input_adps)

        # Magnitudes of the the two sets (by atom)
        output_mags = output_adps[:,0:3].mean(axis=1)
        input_mags = input_adps[:,0:3].mean(axis=1)

        # If output < input, leave as is.
        leave_sel = (output_mags <= input_mags)
        output_mags[leave_sel] = 1.
        input_mags[leave_sel] = 1.

        # Multipliers needed to scale output -> input
        multipliers = (input_mags / output_mags)

        assert (multipliers <= 1.).all()

        rescaled_adps = multipliers.reshape((len(multipliers),1)) * output_adps
        rescaled_adps = flex.sym_mat3_double(list(map(tuple, rescaled_adps)))

        return rescaled_adps

