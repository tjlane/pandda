


from cctbx.uctbx import unit_cell
from cctbx.sgtbx import space_group

from mmtbx import masks

point = (3,4,5)

uc = unit_cell('3,4,5,90,90,90')
sg = space_group('P1')
grid_spacing = 0.1

radius = 1
am = masks.atom_mask(unit_cell=uc, space_group=sg, gridding_n_real=(30,40,50), solvent_radius=1.8)
