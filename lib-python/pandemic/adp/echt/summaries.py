# Summaries that may be used in the future...

def format_no_fail(f_str, val):
    try:
        s = ('{'+f_str+'}').format(val)
    except:
        s = str(val)
    return s

def input_data_summary(
        model_object,
        isotropic_mask = None,
        ):

    mo = model_object

    if isotropic_mask is None:
        disorder_model = 'mixed'
        n_iso = isotropic_mask.iselection().size()
        n_ani = isotropic_mask.size() - n_iso
        n_av = (n_iso + 6*n_ani) / float(isotropic_mask.size())
    else:
        disorder_model = 'anisotropic'

    s = ''
    s += '\n> Input datasets:'
    s += '\nNumber of datasets: {}'.format(mo.n_datasets)
    s += '\nNumber of atoms (in each structure): {}'.format(mo.n_atoms)
    s += '\nInput Disorder Type: {}'.format(disorder_model)
    if disorder_model == 'mixed':
        s += '\n\t{:d} anisotropic atoms ({:5.2%})'.format(n_ani, n_ani/float(mo.n_atoms))
        s += '\n\t{:d} isotropic atoms   ({:5.2%})'.format(n_iso, n_iso/float(mo.n_atoms))
        s += '\nNumber of uij values per atom (average): {:.2f}'.format(n_av)
        s += '\nNumber of uij values (total): {}'.format(n_av*mo.n_atoms)
    else:
        s += '\nNumber of uij values per atom (average): {:.2f}'.format(6)
        s += '\nNumber of uij values (total): {}'.format(6*mo.n_atoms)

    return s

def parameterisation_optimisation_summary():

    n_adp_params = (n_iso + 6*n_ani)
    n_tls_params_real = [sum([o.tls_parameters.n_params(free=True, non_zero=True) for o in gs]) for gs in mo.tls_objects]

    s += '\nNumber of model parameters per group: {} x (20 + {}) = {}'.format(
            mo.n_modes, mo.n_datasets, n_tls_params_per_group,
            )
    s += '\nTotal TLS Parameters: {}'.format(n_tls_params)
    s += '\nNon-zero tls parameters:'
    for i_l, l in mo.tls_level_names:
        s += '\n\tLevel {} ({}): {} groups'.format(i_l+1, l, n_tls_params_real[i_l])
    s += '\nTotal non-zero tls parameters: {}'.format(sum(n_tls_params_real))

    s += '\nNumber of parameters in {} level: 1*{} + 6*{} = {}'.format(
            mo.adp_level_name, n_iso, n_ani, n_adp_params,
            )
    s += '\nNumber of parameters in {} level: {}'.format(mo.adp_level_name, 0)
    s += '\nNumber of model parameters (total): {}'.format(n_tls_params + n_adp_params)

    s += '\n> Optimisation:'
    s += '\nDatasets used for TLS and residual optimisation: {} of {}'.format(sum(self.dataset_mask), len(self.dataset_mask))
    s += '\nAtoms used for TLS optimisation: {} of {}'.format(sum(self.atomic_mask), len(self.atomic_mask))
    s += '\nDataset Weighting: {}'.format(self.params.optimisation.dataset_weights)
    s += '\nAtomic Weighting: {}'.format(self.params.optimisation.atom_weights)
    s += '\nAtomic weights renormalised by dataset: {}'.format("Yes" if self.params.optimisation.renormalise_atom_weights_by_dataset else "No")
    s += '\n'
    s += '\n> Parameterisation:'
    s += '\nNumber of parameters per atom:'
    s += '\n... in input data (average over all atoms): {}'.format(format_no_fail(':2.2f', self.n_input_params_per_atom()))
    s += '\n... in fitted model: {:5.2f}'.format(self.n_params_per_atom(non_zero=False))
    s += '\n... in fitted model (non-zero): {:5.2f}'.format(self.n_params_per_atom(non_zero=True))
    s += '\nNon-zero model parameters / number of input data: {}'.format(format_no_fail(':5.2%', self.parameter_ratio_gain(non_zero=True)))
    s += '\n'
    return s
