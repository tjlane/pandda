from cctbx import maptbx
import iotbx.pdb
from iotbx.reflection_file_utils import extract_miller_array_from_file
from scitbx.math import superpose
from mmtbx.maps.superpose import transform_map_by_lsq_fit

from Giant.Xray.Maps.Utils import get_fft_map_from_f_obs_and_structure

# Reference data
xray_inp1 = iotbx.pdb.input(source_info=None,lines=open('/home/npearce/LIGANDTEST/x529/1-apo/apo-refmac.pdb','r').read())
xray_str1 = xray_inp1.xray_structure_simple()
ref_str = xray_str1
ref_arr = extract_miller_array_from_file('/home/npearce/LIGANDTEST/x529/1-apo/apo-refmac.mtz','F,SIGF')

# Data to be aligned
xray_inp2 = iotbx.pdb.input(source_info=None,lines=open('/home/npearce/LIGANDTEST/x575/apo.pdb','r').read())
xray_str2 = xray_inp2.xray_structure_simple()
new_str = xray_str2
new_arr = extract_miller_array_from_file('/home/npearce/LIGANDTEST/x575/apo.mtz','F,SIGF')

# Extract the map data
ref_map = get_fft_map_from_f_obs_and_structure(ref_arr, xray_str1)
new_map = get_fft_map_from_f_obs_and_structure(new_arr, xray_str2)

# Extract the atoms to generate the alignment matrix
ref_sites = ref_str.sites_cart()
new_sites = new_str.sites_cart()

print 'Size',ref_sites.size()

# Generate rotation and translation for alignment
lsq_fit_obj = superpose.least_squares_fit(ref_sites, new_sites)

print 'R:\t', '\n\t'.join([' '.join(map(str,l)) for l in lsq_fit_obj.r.as_list_of_lists()])
print 'T:\t', '\n\t'.join([' '.join(map(str,l)) for l in lsq_fit_obj.t.as_list_of_lists()])

# Superpose data
new_map_superposed = maptbx.superpose_maps(
    unit_cell_1        = new_map.unit_cell(),
    unit_cell_2        = ref_map.unit_cell(),
    map_data_1         = new_map.real_map_unpadded(),
    n_real_2           = ref_map.n_real(),
    rotation_matrix    = lsq_fit_obj.r.elems,
    translation_vector = lsq_fit_obj.t.elems)






