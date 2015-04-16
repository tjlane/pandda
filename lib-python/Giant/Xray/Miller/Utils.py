from scitbx.array_family import flex

def apply_symmetry_to_miller_array(chg_miller, ref_miller):
    """Applies the crystal symmetry, spacegroup info, and unit cell from ref_miller to change_miller - returns a miller array"""

    new_miller_set = chg_miller.set(crystal_symmetry=ref_miller.crystal_symmetry(), unit_cell=ref_miller.unit_cell(), space_group_info=ref_miller.space_group_info())
    new_miller_array = new_miller_set.array(data=chg_miller.data(), sigmas=chg_miller.sigmas())

    assert new_miller_array.is_similar_symmetry(ref_miller), 'Failure in applying symmetry from reference miller array'

    return new_miller_array

def check_miller_arrays_compatible_for_scaling(miller1, miller2, relative_length_tolerance=0.1, absolute_angle_tolerance=5):
    """Takes miller arrays and checks that they can be scaled together"""

    if not miller1.is_similar_symmetry(miller2, relative_length_tolerance=relative_length_tolerance, absolute_angle_tolerance=absolute_angle_tolerance):
        return False

    return True

def scale_miller_array_to_reference(chg_arr, ref_arr, by_bins=True, binning=None, n_bins=10):
    """Takes miller array, and scales it to another miller array"""

    if by_bins:
        # Apply or setup bins
        if binning:
            ref_binner = ref_arr.use_binning(binning)
            chg_binner = chg_arr.use_binning(binning)
        else:
            ref_binner = ref_arr.setup_binner(n_bins=n_bins)
            chg_binner = chg_arr.use_binning_of(ref_arr)

        # Create a blank copy of the chg_arr to be filled with the scaled values
        scl_arr = chg_arr.array(data=flex.double(chg_arr.size()), sigmas=flex.double(chg_arr.size()))

        # Scale!
        for i in ref_binner.range_used():
            # Select the indices in the bin
            ref_sel_bin = ref_binner.selection(i)
            chg_sel_bin = chg_binner.selection(i)

            # Select the reflections in the bin
            ref_arr_bin = ref_arr.select(ref_sel_bin)
            chg_arr_bin = chg_arr.select(chg_sel_bin)

            # Print raw ratio or means
#            print('<ref bin>/<scl bin> = {!s}'.format(ref_arr_bin.mean()/chg_arr_bin.mean()))

            # Scale the data in the bin - note it is the argument of 'scale' that gets scaled
            chg_arr_bin_scl = ref_arr_bin.scale(chg_arr_bin)

            # Add this scaled data to the new array
            scl_arr.data().set_selected(chg_sel_bin, chg_arr_bin_scl.data())

    else:
        # Simple scaling
        scl_arr = ref_arr.scale(chg_arr)

    scl_binner = scl_arr.use_binning_of(ref_arr)

    return scl_arr

def scale_amplitude_arrays_via_intensities(chg_arr_a, ref_arr_a, n_bins=20):
    """Take two amplitude arrays, converts them to intensities, scale them, then convert back to amplitudes"""

    assert ref_arr_a.is_similar_symmetry(chg_arr_a)

    # Convert both to intensity arrays for scaling
    ref_arr_i = ref_arr_a.as_intensity_array()
    ref_arr_i.set_observation_type_xray_intensity()
    chg_arr_i = chg_arr_a.as_intensity_array()
    chg_arr_i.set_observation_type_xray_intensity()

    # Scale mill2 to mill1
    chg_arr_i_scale = scale_miller_array_to_reference(chg_arr=chg_arr_i, ref_arr=ref_arr_i, n_bins=n_bins)
    chg_arr_i_scale.set_observation_type_xray_intensity()

    # Convert back to amplitudes
    chg_arr_a_scale = chg_arr_i_scale.as_amplitude_array()
    chg_arr_a_scale.set_observation_type_xray_amplitude()

    assert chg_arr_a_scale.is_similar_symmetry(chg_arr_a)

#    print('Correlation between unscaled and scaled amplitudes: {!s}'.format(round(chg_arr_a.correlation(chg_arr_a_scale).coefficient(),3)))

    return chg_arr_a_scale








