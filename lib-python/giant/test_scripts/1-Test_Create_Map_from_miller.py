
from iotbx.reflection_file_utils import extract_miller_array_from_file
import mmtbx.f_model
import iotbx.pdb
from cctbx import maptbx

# From Test_Read_N_Scale
obs_arr = frag_mill_scale

xray_inp = iotbx.pdb.input(source_info=None,lines=open('/home/npearce/LIGANDTEST/x575/apo.pdb','r').read())
xray_str = xray_inp.xray_structure_simple()



# Do things with maps!

fmodel = mmtbx.f_model.manager(f_obs = obs_arr, xray_structure = xray_str)

density_map = fmodel.electron_density_map(update_f_part1=False)

map_coeffs = density_map.map_coefficients('Fobs')

fft_map = map_coeffs.fft_map(symmetry_flags=maptbx.use_space_group_symmetry)

real_map = fft_map.real_map()





mtz_dataset = map_coeffs.as_mtz_dataset(column_root_label="Fobs_scaled")

mtz_object = mtz_dataset.mtz_object()

mtz_object.write(file_name='/home/npearce/LIGANDTEST/x575/apo_scaled.mtz')
