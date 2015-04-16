import os, sys, glob

try:
    set_nomenclature_errors_on_read("ignore")
    set_colour_map_rotation_for_map(0)
    set_colour_map_rotation_on_read_pdb(0)
    set_recentre_on_read_pdb(0)
except:
    pass

def load_z_blob_coordinates(z_peak_file):
    """Load up a list of large Z-blob coordinates to allow them to be viewed"""

    z_map_peaks  = open(z_peak_file, 'r').read().split('\n')

    headers = z_map_peaks.pop(0)

    assert headers == 'dtag, rank, blob_peak, blob_size, x, y, z, refx, refy, refz, pdb, mtz'

    z_map_peaks  = [z.split(',') for z in z_map_peaks]

    ranked_list = []

    for i,x in enumerate(z_map_peaks):

        dtag =            x[0]
        rank =            int(x[1])
        peak_val =        float(x[2])
        peak_size =       float(x[3])
        peak_coords =     map(float,x[4:7])
        ref_peak_coords = map(float,x[7:10])

        assert rank == i+1, 'RANKED LIST IN WRONG ORDER'

        hit_dict = {
                    'dtag' : dtag,
                    'rank' : rank,
                    'peak_val' : peak_val,
                    'peak_size': peak_size,
                    'peak_coords' : peak_coords,
                    'ref_peak_coords' : ref_peak_coords
                    }

        ranked_list.append(hit_dict)

    return ranked_list

def close_all():

    for mol in molecule_number_list():
        close_molecule(mol)

def load_new_set(hit_idx):

    ##################
    global work_dir
    global hit_list
    ##################

    close_all()

    current_hit = hit_list[hit_idx]

    d_tag = current_hit['dtag']

    dataset_dir = os.path.join(work_dir, 'interesting_datasets/Dataset-{!s}'.format(d_tag))
    lig_dir = os.path.join(dataset_dir, 'ligand_files')
    statmap_dir = os.path.join(work_dir, 'statistical_maps')
    ref_dir     = os.path.join(work_dir, 'reference')

    ########################################################

    p = read_pdb(os.path.join(dataset_dir, '{!s}-aligned.pdb'.format(d_tag)))

    s = handle_read_ccp4_map(os.path.join(dataset_dir, '{!s}-observed.ccp4'.format(d_tag)), 0)
    set_last_map_contour_level(1)
    set_map_displayed(s, 1)

    z = handle_read_ccp4_map(os.path.join(dataset_dir, '{!s}-z_map_adjusted_normalised.ccp4'.format(d_tag)), 1)
    set_last_map_contour_level(3)
    set_map_displayed(z, 1)

    set_rotation_centre(*current_hit['ref_peak_coords'])

    ########################################################

    # Stat maps - TODO REMOVE AS SHOULD ALWAYS BE LOADED

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

    print '======================>'
    print '======================>'
    print '======================>'
    print 'INPUT FOR DATASET {!s}'.format(d_tag)
    print current_hit
    print '======================>'
    print '======================>'
    print '======================>'

    return [p, s, z]

def print_results():

    ##################
    global d_results
    ##################

    print '======================>'
    print '======================>'
    print '======================>'
    print d_results
    print '======================>'
    print '======================>'
    print '======================>'

def key_mark_interesting():

    global d_results
    global hit_idx
    global hit_list

    d_results[hit_idx] = hit_list[hit_idx]

    print '======================>'
    print 'DATASET {!s} MARKED AS INTERESTING AT {!s}'.format(hit_list[hit_idx]['dtag'], hit_list[hit_idx][ref_peak_coords])
    print '======================>'

def key_binding_next():

    global hit_idx
    global hit_list

    hit_idx += 1

    if hit_idx == len(hit_list):
        print '======================>'
        print 'REACHED END OF LIST'
        print '======================>'
        hit_idx = len(hit_list) - 1
    else:
        print '======================>'
        print 'MOVING TO NEXT HIT: {!s}'.format(hit_idx+1)
        print '======================>'
        load_new_set(hit_idx)

def key_binding_prev():

    global hit_idx

    hit_idx -= 1

    if hit_idx == -1:
        print '======================>'
        print 'REACHED END OF LIST'
        print '======================>'
        hit_idx = 0
    else:
        print '======================>'
        print 'MOVING TO PREV HIT: {!s}'.format(hit_idx+1)
        print '======================>'
        load_new_set(hit_idx)

try:
    add_key_binding("next blob", "s", lambda: key_binding_next())
    add_key_binding("prev blob", "a", lambda: key_binding_prev())
    add_key_binding("mark_interesting", "v", lambda: key_mark_interesting())

    add_key_binding("print results", "b", lambda: print_results())
except:
    pass

################################

#global work_dir
#global hit_list
#global hit_idx

work_dir = os.getcwd()

# FIND THE NEXT OUTPUT FILE
out_idx = 1
out_template = 'coot_reviews_{:04d}.csv'

results_file = os.path.join(work_dir, out_template.format(out_idx))
while os.path.exists(results_file):
    if out_idx > 9999:  break
    out_idx += 1
    results_file = os.path.join(work_dir, out_template.format(out_idx))

print '======================>'
print 'WRITING OUTPUT TO:'
print '======================>'
print results_file
print '======================>'

hit_list = load_z_blob_coordinates(os.path.join(work_dir, 'analyses', 'blob_site_summaries.csv'))
hit_idx = 0

d_results = {}

load_new_set(hit_idx)

#http://strucbio.biologie.uni-konstanz.de/ccp4wiki/index.php/Ncs_rotamer_differences.py

