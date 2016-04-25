from iotbx.reflection_file_utils import extract_miller_array_from_file
from Giant.Xray.Miller.Utils import apply_symmetry_to_miller_array, check_miller_arrays_compatible_for_scaling, scale_miller_array_to_reference, scale_amplitude_arrays_via_intensities

ref_mill = extract_miller_array_from_file('/home/npearce/LIGANDTEST/x529/1-apo/apo-refmac.mtz','F,SIGF')
frag_mill = extract_miller_array_from_file('/home/npearce/LIGANDTEST/x575/apo.mtz','F,SIGF')

assert check_miller_arrays_compatible_for_scaling(ref_mill, frag_mill)

# Apply the unit cell, SG, cell dims from mill1 to mill2 for processing
frag_mill_sym = apply_symmetry_to_miller_array(chg_miller=frag_mill, ref_miller=ref_mill)

# Scale
frag_mill_sym_scale = scale_amplitude_arrays_via_intensities(chg_arr_a=frag_mill_sym, ref_arr_a=ref_mill, n_bins=20)

# Apply the unit cell, SG, cell dims from mill2 BACK to mill2_scale
frag_mill_scale = apply_symmetry_to_miller_array(chg_miller=frag_mill_sym_scale, ref_miller=frag_mill)






