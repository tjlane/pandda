import cctbx.miller

import numpy

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

def calculate_wilson_b_factor(intensities, inverse_resolution, low_res_cutoff=4.0):

    ln_ints = numpy.log(intensities)
    sqr_res = numpy.power(inverse_resolution, 2)

    # Truncate data at low_res_cutoff
    idx_cut = numpy.where(inverse_resolution > (1.0/low_res_cutoff))[0]
    if len(idx_cut)==0: raise Exception('Data does not extend to res_cutoff')
    ln_ints = ln_ints[idx_cut[0]:]
    sqr_res = sqr_res[idx_cut[0]:]

    # polyfit
    scale, offset = numpy.polyfit(x=sqr_res, y=ln_ints, deg=1)

    return (scale, offset)

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
