import os, sys, glob, copy

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

    if z_map_peaks[0].startswith('dtag'):
        headers = z_map_peaks.pop(0).split(', ')

    z_map_peaks  = [z.split(',') for z in z_map_peaks]

    ranked_list = []

    for i,x in enumerate(z_map_peaks):

        if x == ['']:
            print 'BLANK ROW'
            continue

        dtag            = x[headers.index('dtag')]
        rank            = i+1
        event           = int(x[headers.index('event_idx')])
        try:       site = int(x[headers.index('site_idx')])
        except:    site = 0
        peak_val        = float(x[headers.index('z_peak')])
        peak_size       = float(x[headers.index('cluster_size')])
        ref_peak_coords = map(float,x[headers.index('refx'):headers.index('refz')+1])
        peak_coords     = map(float,x[headers.index('x'):headers.index('z')+1])

        assert rank == i+1, 'RANKED LIST IN WRONG ORDER'

        hit_dict = {
                    'dtag' : dtag,
                    'rank' : rank,
                    'event': event,
                    'site' : site,
                    'peak_val' : peak_val,
                    'peak_size': peak_size,
                    'peak_coords' : peak_coords,
                    'ref_peak_coords' : ref_peak_coords
                    }

        ranked_list.append(hit_dict)

    return ranked_list

def save_all_models():
    """Save the models to the dataset modelled_structures directory"""

    global model_dir
    global p

    # Find the next available index for the fitted structure
    fitted_outputs = sorted(glob.glob(os.path.join(model_dir, 'fitted-v*')))
    if fitted_outputs:
        last_idx = int(fitted_outputs[-1].replace('.pdb','')[-4:])
        new_fitted = 'fitted-v{:04d}.pdb'.format(last_idx+1)
    else:
        new_fitted = 'fitted-v0001.pdb'
    new_fitted = os.path.join(model_dir, new_fitted)
    # Save the fitted model to this file
    write_pdb_file(p, new_fitted)
    print '======================>'
    print '======================>'
    print '======================>'
    print 'SAVED PROTEIN MOL {!s} TO {!s}'.format(p, new_fitted)
    print '======================>'
    print '======================>'
    print '======================>'

    # Save the pdb file to the output
    fitted_link = os.path.join(model_dir, 'fitted-current.pdb')
    if os.path.exists(fitted_link):
        assert os.path.islink(fitted_link), 'FILE IS NOT A LINK'
        # Delete the symbolic link
        os.unlink(fitted_link)
    print 'Linking {!s} -> {!s}'.format(new_fitted, fitted_link)
    # Create new link the most recent file
    os.symlink(os.path.basename(new_fitted), fitted_link)

def close_all():
    for mol in molecule_number_list():
        close_molecule(mol)

def merge_ligand_with_protein():

    global p, l

    try:
        merge_molecules([l], p)
    except Exception as err:
        print err

def move_ligand_here():

    global l

    try:
        move_molecule_to_screen_centre(l)
    except Exception as err:
        print err

def load_new_set(hit_idx, direction='forwards'):

    ##################
    global work_dir
    global model_dir
    global hit_list
    ##################

    current_hit = hit_list[hit_idx]

    d_tag = current_hit['dtag']

    dataset_dir = glob.glob(os.path.join(work_dir, 'interesting_datasets/*{!s}*'.format(d_tag)))
    if len(dataset_dir) != 1:
        dataset_dir = glob.glob(os.path.join(work_dir, 'interesting_datasets/*{!s}'.format(d_tag)))
    if len(dataset_dir) != 1:
        dataset_dir = glob.glob(os.path.join(work_dir, 'interesting_datasets/{!s}'.format(d_tag)))

    assert len(dataset_dir) == 1, 'Non-unique dataset directory found: {!s}'.format(d_tag)

    dataset_dir = dataset_dir[0]

    statmap_dir = os.path.join(work_dir, 'statistical_maps')
    ref_dir     = os.path.join(work_dir, 'reference')

    lig_dir = os.path.join(dataset_dir, 'ligand_files')
    model_dir   = os.path.join(dataset_dir, 'modelled_structures')
    if not os.path.exists(model_dir): os.mkdir(model_dir)

    # The most recent model of the protein in the pandda maps
    fitted_model = os.path.join(model_dir, 'fitted-current.pdb')

    if (direction=='forwards') and ('--models-only' in sys.argv) and (not os.path.exists(fitted_model)):
        key_binding_next()
        return
    elif (direction=='backwards') and ('--models-only' in sys.argv) and (not os.path.exists(fitted_model)):
        key_binding_prev()
        return

    ########################################################

    if os.path.exists(fitted_model):
        p = read_pdb(fitted_model)
    else:
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

    try:
        occ_map = glob.glob(os.path.join(dataset_dir, '{!s}-event_{!s}_occupancy_*_map.ccp4'.format(d_tag, current_hit['event'])))[0]
        o = handle_read_ccp4_map(occ_map, 0)
        set_last_map_contour_level_by_sigma(1)
        set_map_displayed(o, 1)
    except:
        o = z

    set_rotation_centre(*current_hit['ref_peak_coords'])

    ########################################################

    # Stat maps - TODO REMOVE AS SHOULD ALWAYS BE LOADED

    r = read_pdb(os.path.join(ref_dir, 'reference.symmetry.pdb'))
    set_mol_displayed(r, 0)

#    mean = handle_read_ccp4_map(os.path.join(statmap_dir, 'mean_map.ccp4'), 0)
#    set_map_colour(mean, 0.5, 0.5, 0.5)
#    set_last_map_contour_level(1)
#    set_map_displayed(mean, 1)

    ########################################################

    set_scrollable_map(o)
    set_imol_refinement_map(o)

    ########################################################

    # Ligand Files

    lig_files = glob.glob(os.path.join(lig_dir, '*'))
    lig_pdbs = [f for f in lig_files if f.endswith('.pdb')]
    lig_cifs = [f for f in lig_files if f.endswith('.cif')]
    if (len(lig_pdbs) == 1) and (len(lig_cifs) == 1):
        l = handle_read_draw_molecule_and_move_molecule_here(lig_pdbs[0])
        l_dict = read_cif_dictionary(lig_cifs[0])
        if os.path.exists(fitted_model):
            set_mol_displayed(l, 0)
        else:
            set_mol_displayed(l, 1)
    else:
        l = None

    print '======================>'
    print '======================>'
    print '======================>'
    print 'INPUT FOR DATASET {!s}'.format(d_tag)
    print current_hit
    print '======================>'
    print '======================>'
    print '======================>'
    print 'HIT {!s} of {!s}'.format(hit_idx+1, len(hit_list))
    print '======================>'
    print '======================>'
    print '======================>'

    return p, l

def key_binding_modelled_next():

    save_all_models()
    key_binding_next()

def key_binding_next():

    global hit_idx
    global hit_list
    global p, l

    hit_idx += 1

    if hit_idx == len(hit_list):
        print '======================>'
        print 'REACHED END OF LIST'
        print '======================>'
        hit_idx = 0
    else:
        print '======================>'
        print 'MOVING TO NEXT HIT: {!s}'.format(hit_idx+1)
        print '======================>'
        close_all()
        p, l = load_new_set(hit_idx, direction='forwards')

def key_binding_prev():

    global hit_idx
    global p, l

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
        close_all()
        p, l = load_new_set(hit_idx, direction='backwards')

try:
    add_key_binding("next blob", "s", lambda: key_binding_modelled_next())
    add_key_binding("prev blob", "a", lambda: key_binding_prev())
    add_key_binding("skip blob", "v", lambda: key_binding_next())

    add_key_binding("merge ligand with protein", "m", lambda: merge_ligand_with_protein())

    coot_toolbar_button("MERGE", "merge_ligand_with_protein()")
    coot_toolbar_button("MOVE", "move_ligand_here()")

    add_key_binding("print results", "b", lambda: print_results())
except Exception as err:
    print err

################################

#global work_dir
#global hit_list
#global hit_idx
#global p, l


work_dir = os.getcwd()


#######################################################################################
# Process args
cut_idx = [i for i,a in enumerate(sys.argv) if a.endswith('inspect.py')]
assert len(cut_idx) == 1
args = copy.copy(sys.argv)[cut_idx[0]+1:]
#######################################################################################

## FIND THE NEXT OUTPUT FILE
#out_idx = 1
#out_template = 'coot_reviews_{:04d}.csv'
#
#results_file = os.path.join(work_dir, out_template.format(out_idx))
#while os.path.exists(results_file):
#    if out_idx > 9999:  break
#    out_idx += 1
#    results_file = os.path.join(work_dir, out_template.format(out_idx))
#
#print '======================>'
#print 'WRITING OUTPUT TO:'
#print '======================>'
#print results_file
#print '======================>'

#######################################################################################
hit_list = load_z_blob_coordinates(os.path.join(work_dir, 'analyses', 'identified_sites.csv'))
#######################################################################################


#######################################################################################
# Select datasets to view, or all
dataset_selection = [a for a in args if a.startswith('d=')]
if dataset_selection:
    dataset_selection = dataset_selection[0][2:].split(',')
    print 'ONLY SELECTING DATASETS: {!s}'.format(dataset_selection)
    print '{!s} HITS BEFORE FILTERING'.format(len(hit_list))
    hit_list = [h for h in hit_list if h['dtag'] in dataset_selection]
    print '{!s} HITS AFTER FILTERING'.format(len(hit_list))
else:
    print '{!s} HITS LOADED'.format(len(hit_list))
#######################################################################################


#######################################################################################
# Select a startpoint for the hit list
start_point = [a for a in args if a.startswith('start=')]
if start_point:
    start_point = int(start_point[0][6:])
    assert start_point > 0
    assert start_point > 0
    print 'STARTING FROM HIT: {!s}'.format(start_point)
    hit_idx = start_point-1
else:
    hit_idx = 0
#######################################################################################


#######################################################################################

d_results = {}

p,l = load_new_set(hit_idx)

#http://strucbio.biologie.uni-konstanz.de/ccp4wiki/index.php/Ncs_rotamer_differences.py

