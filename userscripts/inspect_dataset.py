import os, sys, glob

try:
    set_nomenclature_errors_on_read("ignore")
    set_colour_map_rotation_for_map(0)
    set_colour_map_rotation_on_read_pdb(0)
    set_recentre_on_read_pdb(0)
except:
    pass

dataset_dir = os.getcwd()

base_dataset_dir = os.path.basename(dataset_dir)

if '-' in base_dataset_dir:
    d_tag = base_dataset_dir.split('-')[-1]
else:
    d_tag = base_dataset_dir

print 'DTAG:', d_tag
assert d_tag, d_tag

ref_dir     = os.path.join(dataset_dir, '../../reference')
statmap_dir = os.path.join(dataset_dir, '../../statistical_maps')
lig_dir     = os.path.join(dataset_dir, 'ligand_files')

########################################################

p = read_pdb(os.path.join(dataset_dir, '{!s}-aligned.pdb'.format(d_tag)))

s = handle_read_ccp4_map(os.path.join(dataset_dir, '{!s}-observed.ccp4'.format(d_tag)), 0)
set_last_map_contour_level(1)
set_map_displayed(s, 0)

z = handle_read_ccp4_map(os.path.join(dataset_dir, '{!s}-z_map_adjusted_normalised.ccp4'.format(d_tag)), 1)
set_last_map_contour_level(3)
set_map_displayed(z, 1)

d = handle_read_ccp4_map(os.path.join(dataset_dir, '{!s}-mean_diff.ccp4'.format(d_tag)), 1)
set_last_map_contour_level_by_sigma(3)
set_map_displayed(d, 0)

occ_maps = glob.glob(os.path.join(dataset_dir, '{!s}-occupancy-*.ccp4'.format(d_tag)))
for o_map in occ_maps:
    o = handle_read_ccp4_map(o_map, 0)
    set_last_map_contour_level_by_sigma(1)
    set_map_displayed(o, 1)

########################################################

r = read_pdb(os.path.join(ref_dir, 'reference.symmetry.pdb'))
set_mol_displayed(r, 0)

mean = handle_read_ccp4_map(os.path.join(statmap_dir, 'mean_map.ccp4'), 0)
set_map_colour(mean, 0.5, 0.5, 0.5)
set_last_map_contour_level(1)
set_map_displayed(mean, 1)

########################################################

set_scrollable_map(z)

########################################################

# Ligand Files

lig_files = glob.glob(os.path.join(lig_dir, '*'))
lig_pdbs = [f for f in lig_files if f.endswith('.pdb')]
lig_cifs = [f for f in lig_files if f.endswith('.cif')]
if (len(lig_pdbs) == 1) and (len(lig_cifs) == 1):
    l = handle_read_draw_molecule_and_move_molecule_here(lig_pdbs[0])
    l_dict = read_cif_dictionary(lig_cifs[0])
    set_mol_displayed(l, 0)

