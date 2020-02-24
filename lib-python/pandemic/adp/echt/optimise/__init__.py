import os, copy, collections, functools
from libtbx import adopt_init_args, group_args
from bamboo.common.logs import Log
from pandemic.adp.utils import show_file_dict
import numpy

def optimisation_summary_from_target_evaluator(target_evaluator):

        te = target_evaluator

        """Print the optimisation values and weights"""
        s = te.log._bar()
        s += '\nOptimisation Summary: {}\n'.format(self.label)
        s += te.log._bar()
        s += '\nOptimisation Least Squares: {:5.3f}'.format(te.optimisation_lst_sqr)
        s += '\nOptimisation Penalty: {:5.3f}'.format(te.optimisation_penalty)
        s += '\nOptimisation Target: {:5.3f}'.format(te.optimisation_lst_sqr+te.optimisation_penalty)
        s += '\n'+self.log._bar()
        s += '\nModels used for optimisation (and total weights):'
        wts = self._op.wgts.mean(axis=1)
        assert len(wts) == len(self.get_dataset_mask())
        for i_rel, i_abs in enumerate(self.get_dataset_mask()):
            s += '\n\t{:>3d} -> {:.3f}'.format(i_abs,wts[i_rel])
        s += '\n'+te.log._bar()
        if show: te.log(s)
        return s


class LevelTargetUijCalculator:


    def __init__(self,
            uij_target,
            ):
        adopt_init_args(self, locals())

    def __call__(self,
            uij_fitted,
            ignore_level=None,
            ):
        arr = uij_fitted.copy()
        # Zero the current level (as should not be included in target function!)
        if ignore_level is not None:
            arr[ignore_level] = 0.0
        # Sum over all other levels to find total currently modelled
        arr_sum = numpy.sum(arr, axis=0)
        assert arr_sum.shape == self.uij_target.shape
        # Subtract currently modelled from input data
        return self.uij_target - arr_sum


class InitialiseEchtLevel:


    def __init__(self,
            starting_t_matrix = (1.,1.,1.,0.,0.,0.),
            starting_l_matrix = None,
            starting_s_matrix = None,
            starting_amplitudes = None,
            ):
        adopt_init_args(self, locals())

    def __call__(self,
            model_object,
            i_level,
            ):
        tls_objects = model_object.tls_objects[i_level]
        for tls_obj in tls_objects:
            for mode in tls_obj.tls_parameters:
                if self.starting_t_matrix is not None:
                    mode.matrices.set(values=self.starting_t_matrix, component_string='T')
                if self.starting_l_matrix is not None:
                    mode.matrices.set(values=self.starting_l_matrix, component_string='L')
                if self.starting_s_matrix is not None:
                    mode.matrices.set(values=self.starting_s_matrix, component_string='S')

    def summary(self):
        s = ''
        if self.starting_t_matrix is not None:
            s += ('> Setting T matrices to {}\n'.format(self.starting_t_matrix))
        if self.starting_l_matrix is not None:
            s += ('> Setting L matrices to {}\n'.format(self.starting_l_matrix))
        if self.starting_s_matrix is not None:
            s += ('> Setting S matrices to {}\n'.format(self.starting_s_matrix))
        return s


class OptimiseEchtModel:


    def __init__(self,
            optimise_tls_function,
            optimise_adp_function,
            optimise_level_amplitudes_function,
            n_cycles = 1,
            n_cpus = 1,
            verbose = False,
            log = None,
            ):
        if log is None: log = Log()
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

        self.initialise_level = InitialiseEchtLevel()

        from pandemic.adp.echt.optimise.repair import RepairEchtModel
        self.repair_model = RepairEchtModel()

        from pandemic.adp.echt.optimise.sanitise import SanitiseEchtModel
        self.sanitise_model = SanitiseEchtModel(
            tls_matrices_eps = optimise_tls_function.eps_values.tls_matrices_eps,
            tls_amplitudes_eps = optimise_tls_function.eps_values.tls_amplitudes_eps,
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

        self.log('\nOptimising TLS models for {} levels'.format(model_object.n_tls_levels) \
            + ' (+ {} level)'.format(model_object.adp_level_name)*(self.optimise_adp_level is not None))

        if uij_isotropic_mask is None:
            uij_isotropic_mask_or_not = lambda x: x
        else:
            uij_isotropic_mask_or_not = uij_isotropic_mask

        # Class for calculating target for each optimisation
        calculate_uij_target = LevelTargetUijCalculator(uij_target)

        # Prefill all of the arguments to the amplitude optimisation function (called from multiple places with same arguments)
        optimise_amplitudes = functools.partial(
            self.optimise_level_amplitudes,
            uij_target = uij_target,
            uij_target_weights = uij_target_weights,
            uij_isotropic_mask = uij_isotropic_mask,
            level_group_tree  = level_group_tree,
            max_recursions = None,
            # dataset mask TODO
        )

        # Initial amplitude optimisation
        self.log.subheading('Macrocycle {}: '.format(tracking_object.n_cycle)+'Optimising inter-level amplitudes')
        self.repair_model(model_object)
        self.log('Optimising groups of amplitudes for all independent hierarchies')
        optimise_amplitudes(model_object=model_object)

        # Record the object before optimisation
        tracking_object.update(
            uij_target = uij_target,
            uij_lvl = uij_isotropic_mask_or_not(model_object.uijs()),
            step = 'start',
            i_level = range(model_object.n_levels),
            write_graphs = True,
            )

        # "micro-cycles"
        for i_sub_cycle in xrange(self.n_cycles):

            # Iterate through the TLS levels of the fitting
            for i_level in xrange(model_object.n_tls_levels):

                if (tracking_object.n_cycle == 1) and (i_sub_cycle == 0):

                    # Set all e.g. T matrices to diag(1,1,1)
                    self.log.subheading('Initialising values in Level {}'.format(i_level+1))
                    self.log(self.initialise_level.summary())
                    self.initialise_level(model_object, i_level=i_level)

                    # Initial amplitude optimisation to get them away from zero
                    self.log('Optimising groups of amplitudes for all independent hierarchies')
                    optimise_amplitudes(model_object=model_object)

                    if (self.verbose is True):
                        # Update tracking
                        tracking_object.update(
                            uij_target = uij_target,
                            uij_lvl = uij_isotropic_mask_or_not(model_object.uijs()),
                            step = '{}-{} (init-l{})'.format(tracking_object.n_cycle, i_sub_cycle+1, i_level+1),
                            i_level = range(model_object.n_levels),
                            write_graphs = False,
                            )

                self.log.subheading('Macrocycle {}-{}: '.format(tracking_object.n_cycle, i_sub_cycle+1)+'Fitting TLS Groups (level {} - {})'.format(i_level+1, model_object.tls_level_names[i_level]))

                # Optimise the groups for one level
                self.repair_model(model_object)
                model_object.tls_objects[i_level] = self.optimise_tls_level(
                    uij_target = calculate_uij_target(uij_fitted=uij_isotropic_mask_or_not(model_object.uijs()), ignore_level=i_level),
                    uij_target_weights = uij_target_weights,
                    uij_isotropic_mask = uij_isotropic_mask,
                    uij_optimisation_mask = uij_optimisation_mask,
                    tls_groups = model_object.tls_objects[i_level],
                    tls_selections = model_object.tls_selections[i_level],
                    )

                # Optimise the amplitudes between levels
                if self.optimise_level_amplitudes is not None:
                    self.log.subheading('Macrocycle {}-{}: '.format(tracking_object.n_cycle, i_sub_cycle+1)+'Optimising inter-level amplitudes')
                    self.repair_model(model_object)
                    self.log('Optimising groups of amplitudes for all independent hierarchies')
                    optimise_amplitudes(model_object=model_object)

                if (self.verbose is True):
                    # Update tracking
                    tracking_object.update(
                        uij_target = uij_target,
                        uij_lvl = uij_isotropic_mask_or_not(model_object.uijs()),
                        step = '{}-{} (l{})'.format(tracking_object.n_cycle, i_sub_cycle+1, i_level+1),
                        i_level = range(model_object.n_levels),
                        write_graphs = False,
                        )

            if self.optimise_adp_level is not None:

                if model_object.uijs()[..., 0:3].mean(axis=-1).max() == 0.0:
                    self.log.subheading('Macrocycle {}-{}: '.format(tracking_object.n_cycle, i_sub_cycle+1)+'Not optimising ADPs as rest of model is zero')
                else:
                    # Fit the atomic level
                    self.log.subheading('Macrocycle {}-{}: '.format(tracking_object.n_cycle, i_sub_cycle+1)+'Optimising atomic Uijs')

                    # Update the target uij by subtracting contributions from other levels
                    self.repair_model(model_object)
                    model_object.adp_values = self.optimise_adp_level(
                        uij_values = model_object.adp_values,
                        uij_target = calculate_uij_target(uij_fitted=model_object.uijs(), ignore_level=i_level),
                        uij_target_weights = uij_target_weights,
                        uij_isotropic_mask = None,
                        uij_optimisation_mask = uij_optimisation_mask,
                        )

                    # Optimise the amplitudes between levels
                    self.log.subheading('Macrocycle {}-{}: '.format(tracking_object.n_cycle, i_sub_cycle+1)+'Optimising inter-level amplitudes')
                    self.repair_model(model_object)
                    self.log('Optimising groups of amplitudes for all independent hierarchies')
                    optimise_amplitudes(model_object=model_object)

            # Update tracking
            tracking_object.update(
                uij_target = uij_target,
                uij_lvl = uij_isotropic_mask_or_not(model_object.uijs()),
                step = '{}-{}'.format(tracking_object.n_cycle, i_sub_cycle+1),
                i_level = range(model_object.n_levels),
                write_graphs = False,
                )

            # Ensure no negative eigenvalues
            self.repair_model(model_object)

            # Reset null groups in the model ready for next cycle
            self.sanitise_model(model_object)

        # Update tracking
        tracking_object.update(
            uij_target = uij_target,
            uij_lvl = uij_isotropic_mask_or_not(model_object.uijs()),
            step = 'end',
            i_level = range(model_object.n_levels),
            write_graphs = True,
            )

        return model_object


class UpdateOptimisationFunction:


    level_amplitude_string = 'level amplitudes weights'

    show_file_dict = show_file_dict

    def __init__(self,
        gradient_optimisation_decay_factor,
        optimisation_weights_to_update,
        model_optimisation_function,
        plotting_object = None,
        verbose = True,
        log = None,
        ):
        if log is None: log = Log()
        start_values = group_args(
            level_amplitudes_optimisation_weights = copy.deepcopy(model_optimisation_function.optimise_level_amplitudes.optimisation_weights),
            )
        history = collections.OrderedDict()
        result = group_args(
            output_files = None,
            )
        adopt_init_args(self, locals(), exclude=('model_optimisation_function',))

    def update(self,
        model_optimisation_function,
        n_cycle,
        ):

        assert n_cycle >= 1

        decay_factor = self.gradient_optimisation_decay_factor
        total_decay_factor = decay_factor ** (n_cycle-1)

        # Update parameters
        self.log.subheading('Updating Level Amplitude Optimisation Weights')
        self.log('> Cycle {}\n'.format(n_cycle))
        self.log('> Total decay factor (relative to starting values): {}\n'.format(total_decay_factor))

        start_weights = self.start_values.level_amplitudes_optimisation_weights
        new_weights_dict = {}
        for k in self.optimisation_weights_to_update:
            sta_v = getattr(start_weights, k)
            new_v = sta_v / total_decay_factor
            self.log('Updating {} = {} -> {}'.format(k, sta_v, new_v))
            new_weights_dict[k] = new_v
        # Create new object, overriding the selected weights
        new_weights = start_weights.transfer_from_other(new_weights_dict)
        # Replace the original weights object with the new weights object
        model_optimisation_function.optimise_level_amplitudes.optimisation_weights = new_weights

        # Update history
        self.log('\nCurrent values:')
        for k, v in new_weights._asdict().items():
            self.log('{} -> {}'.format(k, v))
            self.history \
                .setdefault(self.level_amplitude_string, collections.OrderedDict()) \
                .setdefault(k, []).append((n_cycle, v))

    def write(self, output_directory):

        self.log.subheading('Writing optimisation weights graphs')

        output_files = collections.OrderedDict()

        for s in [self.level_amplitude_string]:
            prefix = os.path.join(output_directory, s.replace(' ', '_'))
            for variable, values in self.history[s].iteritems():
                x_vals, y_vals = zip(*values)
                filename = prefix+'-{variable}.png'.format(variable=variable)
                fig, axis = self.plotting_object.lineplot(
                    x_vals = x_vals,
                    y_vals = y_vals,
                    title = 'Weight "{}" over cycles'.format(variable),
                    x_label = 'cycle',
                    y_label = 'weight',
                    x_ticks = map(int,x_vals),
                    filename = None, # returns fig and axis
                    marker = '*',
                    markersize = 10,
                    markeredgecolor = 'k',
                    linewidth = 3,
                    )
                # Log scale if data has positive values
                if max(y_vals) > 0.0:
                    try:
                        axis.set_yscale('log', base=10.)
                    except:
                        pass
                self.plotting_object.helper.write_and_close_fig(fig=fig, filename=filename)
                output_files.setdefault(s,collections.OrderedDict())[variable] = filename

        self.show_file_dict(output_files)

        self.result.output_files = output_files

        return output_files

    def as_html_summary(self):
        from pandemic.adp.echt.html.optimise import EchtOptimisationParametersHtmlSummary
        return EchtOptimisationParametersHtmlSummary(self)
