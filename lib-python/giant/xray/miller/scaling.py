from giant.xray.miller.utils import apply_symmetry_to_miller_array, scale_amplitude_arrays_via_intensities

def apply_simple_scaling(miller, ref_miller):
    """Applies a simple scaling algorithm to two similar, but with non-identical unit cells, datasets"""
    miller_sym       = apply_symmetry_to_miller_array(chg_miller=miller, ref_miller=ref_miller)
    miller_sym_scale = scale_amplitude_arrays_via_intensities(chg_arr_a=miller_sym, ref_arr_a=ref_miller, n_bins=20)
    miller_scale     = apply_symmetry_to_miller_array(chg_miller=miller_sym_scale, ref_miller=miller)

    return miller_scale

def divide_miller_indices_into_equal_bins(miller, indices_per_bin):
    """Analyse miller indices and split into equally populated bins"""

    pass

def scale_and_truncate_miller_array(miller, ref_miller, low_index_cutoff, high_index_cutoff, index_bins=None):
    """Apply scaling in resolution bins. Miller arrays are then truncated at low and high miller indices. Scaling bins can be predetermined."""

    pass


if __name__ == '__main__':

    import iotbx.mtz

    m = iotbx.mtz.object('./reference.mtz')

#    c = m.extract_complex('FWT', 'PHWT')
    a=m.as_miller_arrays_dict()[('XDScrystal', 'XDSdataset', 'FWT')]
