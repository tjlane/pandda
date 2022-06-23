import giant.logs as lg
logger = lg.getLogger(__name__)

import copy, math
import numpy
from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure
from scitbx.array_family import flex
from scitbx import matrix
from giant.structure.tls import tls_str_to_n_params

_DEBUG = False


class TLSSimplexGenerator(object):


    _I_sqr = (1.,0.,0.,
              0.,1.,0.,
              0.,0.,1.)

    _sym_mat3_deltas_all = [
                [
                    ( 1., 0., 0., 0., 0., 0.),
                    ( 0., 1., 0., 0., 0., 0.),
                    ( 0., 0., 1., 0., 0., 0.),
                    ( 0., 0., 0., 1., 0., 0.),
                    ( 0., 0., 0., 0., 1., 0.),
                    ( 0., 0., 0., 0., 0., 1.),
                    ],
                [
                    (-1.,-1.,-1., 0., 0., 0.),
                    ( 1.,-1., 0., 0., 0., 0.),
                    ( 0., 1.,-1., 0., 0., 0.),
                    ( 0., 0., 0., 1., 0., 0.),
                    ( 0., 0., 0., 0., 1., 0.),
                    ( 0., 0., 0., 0., 0., 1.),
                    ],
                ]

    def __init__(self,
            vibration_delta,
            libration_delta,
            ):
        adopt_init_args(self, locals())

    def __call__(self,
            tls_matrices,
            components = "TLS",
            ):

        self._validate(components)

        # Number of parameters and start values
        start = tls_matrices.get(components, include_szz=False)
        n = tls_str_to_n_params(components, include_szz=False)
        assert n == len(start)

        # Get the decomposed matrices
        d = tls_matrices.decompose()
        # If decomposition is invalid (fixing mode), obviously can't use the decomposition
        if not d.is_valid(): d = None

        # Iterate through in fixed order
        c_deltas = []
        if 'T' in components:
            deltas = self._sym_mat3_deltas(
                    mat=tls_matrices.T,
                    delta=self.vibration_delta,
                    rot=(d.R_ML if d else None),
                    base_set=0,
                    )
            assert deltas.shape == (6, 6)
            c_deltas.append(deltas)
        if 'L' in components:
            deltas = self._sym_mat3_deltas(
                    mat=tls_matrices.L,
                    delta=self.libration_delta,
                    rot=(d.R_ML if d else None),
                    base_set=0,
                    )
            assert deltas.shape == (6, 6)
            c_deltas.append(deltas)
        if 'S' in components:
            deltas = self._mat3_deltas(
                    mat=tls_matrices.S,
                    delta=self.libration_delta,
                    rot=(d.R_ML if d else None),
                    )
            # Need to truncate Szz
            deltas = deltas[0:8,0:8]
            assert deltas.shape == (8, 8)
            c_deltas.append(deltas)

        # Insert the components into the full simplex
        i_total = 0
        simplex_deltas = numpy.zeros((n+1, n))
        for block in c_deltas:
            n_c = len(block)
            for i, dels in enumerate(block):
                simplex_deltas[1+i_total+i, i_total:i_total+n_c] = dels
            i_total += n_c
        assert i_total == n

        start_simplex = numpy.repeat([start], len(start)+1, axis=0)
        assert start_simplex.shape == simplex_deltas.shape
        start_simplex += simplex_deltas
        return start_simplex

    def _validate(self, components):
        assert set('TLS').intersection(components), 'components must contain T, L or S'

    def _sym_mat3_deltas(self, mat, delta, rot=None, base_set=0):
        """Get deltas caused by changes to a symmetric matrix in a rotated frame"""

        # Extract types of deltas
        assert base_set < len(self._sym_mat3_deltas_all)
        mat_deltas = self._sym_mat3_deltas_all[base_set]

        # Calculate transformation to reference frame
        if rot is None:
            rot = self._I_sqr
        # Create rotation matrix: R maps to diagonalised frame
        R = matrix.sqr(rot)
        R_t = R.transpose()
        # Transform input matrix to the rotated frame
        mat_R = R * matrix.sym(sym_mat3=tuple(mat)) * R_t

        # Output deltas
        n = 6
        deltas = numpy.zeros((n,n))

        for i in range(n):
            # Form the delta matrix in the rotated frame
            d_R = list(mat_deltas[i])
            for i_el in range(n):
                # Get the sign of the current value of the element
                sgn = (1., -1.)[mat_R[i_el] > 0.0]
                d_R[i_el] *= (sgn * delta)
            # Transform back
            d = (R_t * matrix.sym(sym_mat3=tuple(d_R)) * R).as_sym_mat3()
            # Append to output
            deltas[i] = d
        # remark: This check is necessary because providing a numpy float can cause
        # remark: the matrix multiplication to become a numpy array of matrices
        assert deltas.shape == (n,n)
        return deltas

    def _mat3_deltas(self, mat, delta, rot=None):
        """Get deltas caused by changes to a non-symmetric matrix in a rotated frame"""

        # Calculate transformation to reference frame
        if rot is None:
            rot = self._I_sqr
        # Create rotation matrix
        R = matrix.sqr(rot)
        R_t = R.transpose()
        # Transform input matrix to the rotated frame
        mat_R = R * matrix.sqr(tuple(mat)) * R_t

        # Output deltas
        n = 9
        deltas = numpy.zeros((n,n))

        for i in range(n):
            # Form the delta matrix in the rotated frame
            d_R = [0.0,]*n
            # Get the sign of the current value of the element
            sgn = (1., -1.)[mat_R[i] > 0.0]
            d_R[i] = (sgn * delta)
            # Transform back
            d = (R_t * matrix.sqr(tuple(d_R)) * R).as_mat3()
            # Append to output
            deltas[i] = d
        # remark: This check is necessary because providing a numpy float can cause
        # remark: the matrix multiplication to become a numpy array of matrices
        assert (deltas.shape == (n,n))
        return deltas

    def summary(self):
        s = 'Simplex optimisation deltas'
        s += '\nVibration Delta: {}'.format(self.vibration_delta)
        s += '\nLibration Delta: {}'.format(self.libration_delta)
        return s


class OptimiseTLSGroup_Matrices_AdoptValues(object):


    def __init__(self,
            components,
            i_mode,
            ):
        n_parameters = tls_str_to_n_params(components, include_szz=False)
        adopt_init_args(self, locals())

    def __call__(self,
            values,
            tls_parameters,
            include_szz = False,
            ):
        assert len(values) == self.n_parameters
        mde = tls_parameters[self.i_mode]
        mde.matrices.set(
                values = values,
                component_string = self.components,
                include_szz = include_szz,
                )


from pandemic.adp.echt.validate.tls import ValidateTLSMode
class OptimiseTLSGroup_Matrices_TargetEvaluator(object):


    start_target = 1e6

    def __init__(self,
            uij_target,
            target_function,
            other_target_functions,
            multi_dataset_tls_group,
            optimise_components,
            optimise_i_mode,
            uij_isotropic_mask = None,
            ):

        # Init values
        n_call = 0
        best_target = self.start_target
        best_target_terms = None

        # Function for inserting matrix values into multi_dataset_tls_group
        adopt_values = OptimiseTLSGroup_Matrices_AdoptValues(
                components = optimise_components,
                i_mode = optimise_i_mode,
                )

        # Function to check e.g. validity of TLS matrices
        validate = ValidateTLSMode()

        adopt_init_args(self, locals())

    def target(self,
            parameters,
            ):
        self.n_call += 1
        # Insert new values into the tls_parameter object
        self.adopt_values(
                values = parameters,
                tls_parameters = self.multi_dataset_tls_group.tls_parameters,
                )
        if _DEBUG and (self.n_call%20==1): self.print_header_line(parameters)
        # Validate the new combination of parameters and return poor score if invalid
        if not self.validate(tls_mode=self.multi_dataset_tls_group.tls_parameters[self.optimise_i_mode]):
            if _DEBUG: self.print_invalid_line(parameters)
            return 1.1 * self.best_target
        # Extract new uijs
        uij_fitted = self.multi_dataset_tls_group.tls_parameters.uijs(
                sites_carts = self.multi_dataset_tls_group.coordinates,
                origins = self.multi_dataset_tls_group.origins,
                )
        # Make selected atoms isotropic
        if self.uij_isotropic_mask is not None:
            uij_fitted = self.uij_isotropic_mask(uij_fitted)
        # Calculate difference between fitted and target
        uij_deltas = uij_fitted - self.uij_target
        # Calculate target function
        lsq = self.target_function(uij_deltas=uij_deltas)
        # Calculate full target and store if minimum
        target_terms = [lsq] + [f(lsq=lsq, uij_deltas=uij_deltas) for f in self.other_target_functions]
        target = sum(target_terms)
        if target < self.best_target:
            self.best_target = target
            self.best_target_terms = target_terms
        if _DEBUG: self.print_current_line(parameters, target)
        return target

    def print_header_line(self, values):
        header = '[{}] -> ({:^12}, {:^10})'.format(', '.join(['{:>10}'.format('param') for r in values]), 'target', '(best)')
        line = '-'*len(header)
        logger(line)
        logger(header)
        logger(line)

    def print_current_line(self, values, target):
        logger('[{}] -> ({:12.9f}, {:10.7f})'.format(', '.join(['{:+10.7f}'.format(r) for r in values]), target, self.best_target))

    def print_invalid_line(self, values):
        logger('[{}] -> ({:>12}, {:>10})'.format(', '.join(['{:+10.7f}'.format(r) for r in values]), 'UNPHYSICAL', ''))


class OptimiseTLSGroup_Matrices(object):


    evaluator_class = OptimiseTLSGroup_Matrices_TargetEvaluator

    def __init__(self,
            uij_target,
            uij_target_weights,
            simplex_generator,
            target_function,
            other_target_functions,
            convergence_tolerance,
            uij_isotropic_mask = None,
            ):

        # Convert to flex
        uij_target = flex.sym_mat3_double(uij_target.reshape((uij_target.shape[0]*uij_target.shape[1],6)))
        uij_target_weights = flex.double(uij_target_weights.reshape((uij_target_weights.size,)))

        # Copy target functions to avoid interferring with other classes
        target_function = copy.deepcopy(target_function)
        other_target_functions = copy.deepcopy(other_target_functions)

        adopt_init_args(self, locals())

        for f in [self.target_function]+self.other_target_functions:
            f.set_uij_weights(uij_target_weights)

    def __call__(self,
            multi_dataset_tls_group,
            optimise_components,
            optimise_i_mode,
            ):

        # Initialise new instance for each call
        evaluator = self.evaluator_class(
                uij_target = self.uij_target,
                target_function = self.target_function,
                other_target_functions = self.other_target_functions,
                multi_dataset_tls_group = multi_dataset_tls_group,
                optimise_components = optimise_components,
                optimise_i_mode = optimise_i_mode,
                uij_isotropic_mask = self.uij_isotropic_mask,
                )

        # Generate simplex
        mode = multi_dataset_tls_group.tls_parameters[optimise_i_mode]
        starting_simplex = self.simplex_generator(tls_matrices=mode.matrices, components=optimise_components)

        # Optimise these parameters
        from scitbx import simplex
        optimised = simplex.simplex_opt(dimension = len(starting_simplex[0]),
                                        matrix    = list(map(flex.double, starting_simplex)),
                                        evaluator = evaluator,
                                        tolerance = self.convergence_tolerance)

        # Adopt the solution!
        evaluator.adopt_values(
                values = optimised.get_solution(),
                tls_parameters = multi_dataset_tls_group.tls_parameters,
                )

        # Create results object
        results = group_args(
            n_iter = evaluator.n_call,
            target_function_value = evaluator.best_target,
            target_function_values = evaluator.best_target_terms,
            )

        # print evaluator.best_target

        return (multi_dataset_tls_group, results)


