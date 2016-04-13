import cctbx.miller

def extract_structure_factors(mtz_object, ampl_label, phas_label):
    # Get the crystal symmetry from the amplitudes' crystal
    try:
        ampl_col = mtz_object.get_column(ampl_label)
    except:
        raise Exception('Amplitude column not found: {}. Have you specified the right columns?'.format(ampl_label))
    crystal_symmetry = ampl_col.mtz_crystal().crystal_symmetry()
    # Create miller set
    mill_idx = mtz_object.extract_miller_indices()
    mill_set = cctbx.miller.set(crystal_symmetry=crystal_symmetry, indices=mill_idx)
    # Extract amplitudes and phases
    try:
        ampl_com = mtz_object.extract_complex(column_label_ampl=ampl_label, column_label_phi=phas_label)
    except:
        raise Exception('Could not extract structure factors - Amplitudes:{}, Phases:{}. Have you specified the right columns?'.format(ampl_label, phas_label))
    # Convert to miller array
    mill_sfs = mill_set.array(ampl_com.data)
    mill_sfs.set_observation_type_xray_amplitude()
    # Check it's complex
    assert mill_sfs.is_complex_array(), 'STRUCTURE FACTORS SHOULD BE COMPLEX?!'
    # Make non-anomalous
    mill_sfs = mill_sfs.as_non_anomalous_array()
    return mill_sfs
