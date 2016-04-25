from cctbx import maptbx
import iotbx.pdb
from iotbx.reflection_file_utils import extract_miller_array_from_file
from scitbx.math import superpose

from Giant.Xray.Maps.Utils import get_fft_map_from_f_obs_and_structure

ref_arr = extract_miller_array_from_file('/home/npearce/LIGANDTEST/x529/1-apo/apo-refmac.mtz','F,SIGF')

xray_inp1 = iotbx.pdb.input(source_info=None,lines=open('/home/npearce/LIGANDTEST/x529/1-apo/apo-refmac.pdb','r').read())
xray_str1 = xray_inp1.xray_structure_simple()

new_arr = extract_miller_array_from_file('/home/npearce/LIGANDTEST/x575/apo.mtz','F,SIGF')

xray_inp2 = iotbx.pdb.input(source_info=None,lines=open('/home/npearce/LIGANDTEST/x575/apo.pdb','r').read())
xray_str2 = xray_inp2.xray_structure_simple()

ref_map = get_fft_map_from_f_obs_and_structure(ref_arr, xray_str1)
new_map = get_fft_map_from_f_obs_and_structure(new_arr, xray_str2)

search_params = maptbx.peak_search_parameters(max_peaks=30)

ref_peaks = ref_map.peak_search(search_params)
new_peaks = new_map.peak_search(search_params)

ref_sites = ref_peaks.sites()
new_sites = new_peaks.sites()




