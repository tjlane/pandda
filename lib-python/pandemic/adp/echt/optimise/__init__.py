import giant.logs as lg
logger = lg.getLogger(__name__)

import os, copy, functools
from libtbx import adopt_init_args

from .convergence import EchtModelAmplitudeConvergenceChecker

# def optimisation_summary_from_target_evaluator(target_evaluator):

#     te = target_evaluator

#     """Print the optimisation values and weights"""
#     s += '\nOptimisation Summary: {}\n'.format(self.label)
#     s += '\nOptimisation Least Squares: {:5.3f}'.format(te.optimisation_lst_sqr)
#     s += '\nOptimisation Penalty: {:5.3f}'.format(te.optimisation_penalty)
#     s += '\nOptimisation Target: {:5.3f}'.format(te.optimisation_lst_sqr+te.optimisation_penalty)
#     s += '\nModels used for optimisation (and total weights):'
#     wts = self._op.wgts.mean(axis=1)
#     assert len(wts) == len(self.get_dataset_mask())
#     for i_rel, i_abs in enumerate(self.get_dataset_mask()):
#         s += '\n\t{:>3d} -> {:.3f}'.format(i_abs,wts[i_rel])
#     return s


class LevelTargetUijCalculator(object):

    def __init__(self,
            uij_target,
            ):
        adopt_init_args(self, locals())

    def __call__(self,
            uij_fitted,
            ignore_level=None,
            ):
        import numpy
        arr = uij_fitted.copy()
        # Zero the current level (as should not be included in target function!)
        if ignore_level is not None:
            arr[ignore_level] = 0.0
        # Sum over all other levels to find total currently modelled
        arr_sum = numpy.sum(arr, axis=0)
        assert arr_sum.shape == self.uij_target.shape
        # Subtract currently modelled from input data
        return self.uij_target - arr_sum


class OptimiseEchtModel(object):

    debug = False
    
    def __init__(self,
            optimise_tls_function,
            optimise_adp_function,
            optimise_level_amplitudes_function,
            max_n_cycles = 1,
            n_cpus = 1,
            ):
        adopt_init_args(self, locals())

        from pandemic.adp.echt.optimise.tls import OptimiseTLSLevel
        self.optimise_tls_level = OptimiseTLSLevel(
            optimisation_function = self.optimise_tls_function,
            n_cpus = n_cpus,
            )

        if self.optimise_adp_function is not None:
            from pandemic.adp.echt.optimise.uij import OptimiseUijLevel
            self.optimise_adp_level = OptimiseUijLevel(
                optimisation_function = self.optimise_adp_function,
                n_cpus = n_cpus,
                )
        else:
            self.optimise_adp_level = None

        self.optimise_level_amplitudes = optimise_level_amplitudes_function

        # Initialise sanitise
        sanitise_tls_parameters = sanitise_uij_parameters = None
        if optimise_tls_function is not None:
            sanitise_tls_parameters = dict(
                matrix_eps = optimise_tls_function.eps_values.matrix_eps,
                amplitude_eps = optimise_tls_function.eps_values.amplitude_eps,
            )
        if optimise_adp_function is not None:
            sanitise_uij_parameters = dict(
                eps = optimise_adp_function.uij_eps,
                tolerance = optimise_adp_function.uij_tolerance,
            )

        from pandemic.adp.echt.sanitise import SanitiseEchtModel
        self.sanitise_model = SanitiseEchtModel(
            tls_parameters_dict = sanitise_tls_parameters,
            uij_parameters_dict = sanitise_uij_parameters,
            )

    def __call__(self,
            model_object,
            level_group_tree,
            uij_target,
            uij_target_weights,
            uij_isotropic_mask,
            tracking_object,
            uij_optimisation_mask = None,
            ):

        logger('\nOptimising TLS models for {} levels'.format(model_object.n_tls_levels) \
            + ' (+ {} level)'.format(model_object.adp_level_name)*(self.optimise_adp_level is not None))

        if uij_isotropic_mask is None:
            uij_isotropic_mask_or_not = lambda x: x
        else:
            uij_isotropic_mask_or_not = uij_isotropic_mask

        # Class for calculating target for each optimisation
        calculate_uij_target = LevelTargetUijCalculator(uij_target)

        # Prefill all of the arguments to the amplitude optimisation function (called from multiple places with same arguments)
        #
        # All recursions -- and ADPs
        #
        optimise_amplitudes_full = functools.partial(
            self.optimise_level_amplitudes,
            uij_target = uij_target,
            uij_target_weights = uij_target_weights,
            uij_isotropic_mask = uij_isotropic_mask,
            level_group_tree = level_group_tree,
            optimise_atomic_adp_amplitudes = (self.optimise_adp_level is not None),
            max_recursions = None,
            # dataset mask TODO
        )
        #
        # Iterative recursions -- TLS only
        #
        #optimise_amplitudes_iter = functools.partial(
        #    self.optimise_level_amplitudes,
        #    uij_target = uij_target,
        #    uij_target_weights = uij_target_weights,
        #    uij_isotropic_mask = uij_isotropic_mask,
        #    level_group_tree = level_group_tree,
        #    optimise_atomic_adp_amplitudes = (self.optimise_adp_level is not None),
        #    max_recursions = 1,
        #    recursion_direction = 'descending',
        #    # dataset mask TODO
        #)

        # Initial amplitude optimisation
        logger.subheading(
            'Macrocycle {}: '.format(tracking_object.n_cycle)+'Optimising inter-level amplitudes'
            )
        self.sanitise_model(model_object)
        #optimise_amplitudes_iter(model_object=model_object)
        optimise_amplitudes_full(
            model_object = model_object,
            )

        # Record the object before optimisation
        tracking_object.update(
            uijs = uij_isotropic_mask_or_not(model_object.uijs()),
            step = 'start',
            i_level = list(range(model_object.n_levels)),
            )

        is_converged = EchtModelAmplitudeConvergenceChecker()

        # "micro-cycles"
        for i_sub_cycle in range(self.max_n_cycles):

            # Break loop if model is zero
            if self.max_u_iso(model_object) == 0.0:
                logger.subheading(
                    (
                        'Macrocycle {n_cycle}: '
                        'All model amplitudes are zero. Not doing any optimisation...'
                        ).format(
                        n_cycle = tracking_object.n_cycle,
                        )
                    )
                break

            # Iterate through the TLS levels of the fitting
            for i_level in range(model_object.n_tls_levels):

                logger.subheading(
                    (
                        'Macrocycle {n_cycle}-{n_sub_cycle}: '
                        'Fitting TLS Groups (level {level_num} - {level_name})'
                        ).format(
                        n_cycle = tracking_object.n_cycle, 
                        n_sub_cycle = i_sub_cycle+1,
                        level_num = i_level+1,
                        level_name = model_object.tls_level_names[i_level],
                        )
                    )

                model_object.tls_objects[i_level] = self.optimise_tls_level(
                    uij_target = calculate_uij_target(
                        uij_fitted = uij_isotropic_mask_or_not(model_object.uijs()), 
                        ignore_level = i_level,
                        ),
                    uij_target_weights = uij_target_weights,
                    uij_isotropic_mask = uij_isotropic_mask,
                    uij_optimisation_mask = uij_optimisation_mask,
                    tls_groups = model_object.tls_objects[i_level],
                    tls_selections = model_object.tls_selections[i_level],
                    )

                # Optimise the amplitudes between levels
                logger.subheading(
                    (
                        'Macrocycle {n_cycle}-{n_sub_cycle}: '
                        'Optimising inter-level amplitudes'
                        ).format(
                        n_cycle = tracking_object.n_cycle, 
                        n_sub_cycle = i_sub_cycle+1,
                        )
                    )

                self.sanitise_model(model_object)

                optimise_amplitudes_full(
                    model_object = model_object,
                    )

                if (self.debug is True):
                    
                    tracking_object.update(
                        uijs = uij_isotropic_mask_or_not(model_object.uijs()),
                        step = '{n_cycle}-{n_sub_cycle} (l{n_level})'.format(
                            n_cycle = tracking_object.n_cycle, 
                            n_sub_cycle = i_sub_cycle+1, 
                            n_level = i_level+1,
                            ),
                        i_level = list(range(model_object.n_levels)),
                        )

            if self.optimise_adp_level is not None:

                logger.subheading(
                    (
                        'Macrocycle {n_cycle}-{n_sub_cycle}: '
                        'Optimising atomic uijs'
                        ).format(
                        n_cycle = tracking_object.n_cycle, 
                        n_sub_cycle = i_sub_cycle+1,
                        )
                    )

                model_object.adp_values = self.optimise_adp_level(
                    uij_values = model_object.adp_values,
                    uij_target = calculate_uij_target(
                        uij_fitted = uij_isotropic_mask_or_not(model_object.uijs()), 
                        ignore_level = -1,
                        ), # !!! Should this not use the isotropic mask?
                    uij_target_weights = uij_target_weights,
                    uij_isotropic_mask = uij_isotropic_mask, # !!! Should this use the isotropic mask?
                    uij_optimisation_mask = uij_optimisation_mask,
                    )

                # Optimise the amplitudes between levels
                logger.subheading(
                    (
                        'Macrocycle {n_cycle}-{n_sub_cycle}: '
                        'Optimising inter-level amplitudes'
                        ).format(
                        n_cycle = tracking_object.n_cycle, 
                        n_sub_cycle = i_sub_cycle+1,
                        )
                    )

                self.sanitise_model(model_object)

                optimise_amplitudes_full(
                    model_object = model_object,
                    )

            # Update tracking
            tracking_object.update(
                uijs = uij_isotropic_mask_or_not(model_object.uijs()),
                step = '{n_cycle}-{n_sub_cycle}'.format(
                    n_cycle = tracking_object.n_cycle, 
                    n_sub_cycle = i_sub_cycle+1,
                    ),
                i_level = list(range(model_object.n_levels)),
                )

            if is_converged(model_object):
                break

        # Close processes to clear memory, and to make KeyboardInterrupts during this easier to handle (processes are automatically reopened)
        self.close_processes()

        # Optimise the amplitudes again just to make sure...
        logger.subheading(
            (
                'Macrocycle {n_cycle}: '
                'Optimising inter-level amplitudes'
                ).format(
                n_cycle = tracking_object.n_cycle,
                )
            )

        self.sanitise_model(model_object)

        optimise_amplitudes_full(
            model_object = model_object,
            )

        # Update tracking
        tracking_object.update(
            uijs = uij_isotropic_mask_or_not(model_object.uijs()),
            step = 'end',
            i_level = list(range(model_object.n_levels)),
            )

        return model_object

    def max_u_iso(self, model_object):
        return model_object.uijs()[..., 0:3].mean(axis=-1).max()

    def close_processes(self):
        """Terminate processes after optimisation is complete"""
        if self.optimise_tls_level is not None:
            self.optimise_tls_level.run_parallel.close_processes()
        if self.optimise_adp_level is not None:
            self.optimise_adp_level.run_parallel.close_processes()

