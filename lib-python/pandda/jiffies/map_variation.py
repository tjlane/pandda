
import os, sys, glob, copy, time, itertools
import gtk

#=========================================================================
# GTK FUNCTIONS
#=========================================================================

def catchup(block=False):
    while gtk.events_pending():
        gtk.main_iteration(block)

def nonmodal_msg(msg):
    """Display an error window - non-modal"""
    d = gtk.MessageDialog(  type    = gtk.MESSAGE_INFO,
                            message_format = msg )
    d.set_position(gtk.WIN_POS_CENTER)
    d.set_keep_above(True)
    d.show_all()
    catchup()
    return d

def modal_msg(msg):
    """Display an error window - model"""
    d = gtk.MessageDialog(  type    = gtk.MESSAGE_INFO,
                            buttons = gtk.BUTTONS_CLOSE,
                            message_format = msg )
    d.set_position(gtk.WIN_POS_CENTER)
    d.set_keep_above(True)
    d.run()
    d.destroy()

#=========================================================================
# COOT FUNCTIONS
#=========================================================================

def set_main_coot_molecule(i):
    # Have to set colour manually!
    set_molecule_bonds_colour_map_rotation(i, MOL_COLOUR); graphics_draw()
    # Other settings
    set_pointer_atom_molecule(i)
    set_go_to_atom_molecule(i)
    update_go_to_atom_window_on_new_mol()

def coot_customisation():
    set_nomenclature_errors_on_read("ignore")
    set_recentre_on_read_pdb(1)
    set_show_symmetry_master(1)
    set_symmetry_shift_search_size(2)
    set_colour_map_rotation_for_map(0.0)
    set_colour_map_rotation_on_read_pdb(0.0)
    set_colour_map_rotation_on_read_pdb_flag(1)
    set_colour_map_rotation_on_read_pdb_c_only_flag(1)

    add_key_binding("Add ligand",  "a", lambda: solvent_ligands_gui())
    add_key_binding("Add water",   "w", lambda : place_typed_atom_at_pointer("Water"))

def post_coot_windows():
    post_display_control_window()
    post_go_to_atom_window()
    try:
        # Future-available function
        post_delete_item_dialog()
    except:
        pass

#=========================================================================

if __name__=='__main__':

    #############################################################################################
    #
    # CHANGE COOT SETTINGS
    #
    #############################################################################################

    try:
        coot_customisation()
    except:
        pass

    work_dir = os.getcwd()

    reference_pdb = os.path.join(work_dir, 'reference', 'reference.shifted.pdb')
    read_pdb(reference_pdb)

    stat_dir = os.path.join(work_dir, 'reference', 'statistical_maps')
    mean_maps = sorted(glob.glob(os.path.join(stat_dir, '*mean*.ccp4')))

    m1 = mean_maps[0]
    m1_i = handle_read_ccp4_map(m1, 0)
    set_last_map_contour_level(1)
    set_map_displayed(m1_i, 1)

    s1 = m1.replace('mean','sadj')
    s1_i = handle_read_ccp4_map(s1, 0)
    set_last_map_contour_level(0.25)
    set_map_displayed(s1_i, 1)

    b1 = m1.replace('mean','bimo')
    b1_i = handle_read_ccp4_map(b1, 0)
    set_last_map_contour_level(3.0)
    set_map_displayed(b1_i, 1)


