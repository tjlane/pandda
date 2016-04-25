import iotbx.map_tools
from Giant.Xray.Maps.Utils import generate_p1_box_from_size


# Generate p1 box for unit cell
sizes_of_p1_box = ref_unit_cell.orthogonalize((1,1,1))
p1_box = generate_p1_box_from_size(sizes_of_p1_box)
p1_unit_cell = p1_box.unit_cell()

# sample_grid_cart
# new_p1_sampling_unit_cell
# sampled_map_data
# grid_size
# output_filename

iotbx.map_tools.write_ccp4_map(XXX_grid_sites, p1_unit_cell, XXX_map_data, XXX_n_real, XXX_file_name, buffer=0)
