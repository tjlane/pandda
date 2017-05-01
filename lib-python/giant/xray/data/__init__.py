import cctbx.miller

from scitbx.array_family import flex

from giant.stats.optimisation import LinearScaling

def extract_structure_factors(mtz_object, ampl_label, phas_label):
    # Get the crystal symmetry from the amplitudes' crystal
    try:
        ampl_col = mtz_object.get_column(ampl_label)
    except:
        raise Exception('Amplitude column not found: {}. Have you specified the right columns?'.format(ampl_label))
    # Get the symmetry associated with the column
    crystal_symmetry = ampl_col.mtz_crystal().crystal_symmetry()
    # Extract amplitudes and phases
    try:
        sf_com = mtz_object.extract_complex(column_label_ampl=ampl_label, column_label_phi=phas_label)
    except:
        raise Exception('Could not extract structure factors - Amplitudes:{}, Phases:{}. Have you specified the right columns?'.format(ampl_label, phas_label))
    # Convert to miller array
    mill_set = cctbx.miller.set(crystal_symmetry=crystal_symmetry, indices=sf_com.indices)
    mill_sfs = mill_set.array(sf_com.data)
#    mill_sfs.set_observation_type_xray_amplitude()
    mill_sfs = mill_sfs.as_non_anomalous_array()
    assert mill_sfs.is_complex_array(), 'STRUCTURE FACTORS SHOULD BE COMPLEX?!'
    return mill_sfs

def estimate_wilson_b_factor(miller_array, low_res_cutoff=4.0):

    miller_array = miller_array.resolution_filter(d_max=low_res_cutoff).as_intensity_array()

    # Setup binner and extract radial averages
    binner = miller_array.setup_binner(auto_binning=True)
    binned = miller_array.wilson_plot(use_binning=True)
    # Convert to scale
    y_values = flex.log(flex.double(binned.data[1:-1]))
    x_values = flex.pow2(binner.bin_centers(1))
    # Perform scaling
    scl = LinearScaling(x_values   = x_values,
                        ref_values = y_values)

    return -0.5*scl.optimised_values[1]

#def extract_structure_factors(mtz_object, ampl_label, phas_label):
#
#    # Extract matching miller arrays
#    match_arrs = [a for a in mtz_object.as_miller_arrays() if a.info().labels==[ampl_label, phas_label]]
#    if not match_arrs: raise Exception('Could not extract structure factors - Amplitudes:{}, Phases:{}. Have you specified the right columns?'.format(ampl_label, phas_label))
#    assert len(match_arrs) == 1
#    mill_arr = match_arrs[0]
#    assert mill_arr.is_complex_array(), 'STRUCTURE FACTORS SHOULD BE COMPLEX?!'
#    mill_arr = mill_arr.as_non_anomalous_array()
#    return mill_arr
