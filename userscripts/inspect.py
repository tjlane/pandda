
import os, sys, glob, copy
import gtk
import pandas

try:
    set_nomenclature_errors_on_read("ignore")
    set_colour_map_rotation_for_map(0)
    set_colour_map_rotation_on_read_pdb(0)
    set_recentre_on_read_pdb(0)
except:
    pass

class PanddaEvent(object):
    def __init__(self, rank, info, top_dir):

        # Key for the pandas table
        self.index = info.name
        # Dataset tag
        self.dtag = info.name[0]
        # Event number for the dataset
        self.event_idx = int(info.name[1])
        # Position in the ranked list (1 -> n)
        self.rank = rank
        # Site Number (1 -> m)
        self.site_idx  = int(info['site_idx'])
        # Z statistics
        self.z_peak = info['z_peak']
        self.z_mean = info['z_mean']
        self.cluster_size = int(info['cluster_size'])
        # Coordinate information
        self.ref_coords = info[['refx','refy','refz']]
        self.nat_coords = info[['x','y','z']]
        # Find the files for loading
        self.find_file_paths(top_dir=top_dir)

    def find_file_paths(self, top_dir):
        dtag = self.dtag
        # Top directory fo the pandda processing folder
        self.top_dir = top_dir
        # Identify other common directories
        self.stat_map_dir = os.path.join(top_dir, 'statistical_maps')
        self.ref_dir      = os.path.join(top_dir, 'reference')
        # Symmetry contacts
        self.reference_symmetry = os.path.join(self.ref_dir, 'reference.symmetry.pdb')
        # Find the directory for the dataset
        self.dataset_dir = glob.glob(os.path.join(top_dir, 'interesting_datasets/*{!s}*'.format(dtag)))
        if len(self.dataset_dir) != 1:
            self.dataset_dir = glob.glob(os.path.join(top_dir, 'interesting_datasets/*{!s}'.format(dtag)))
        if len(self.dataset_dir) != 1:
            self.dataset_dir = glob.glob(os.path.join(top_dir, 'interesting_datasets/{!s}'.format(dtag)))
        assert len(self.dataset_dir) == 1, 'Non-unique dataset directory found: {!s}'.format(dtag)
        self.dataset_dir = self.dataset_dir[0]

        # Identify dataset subdirectories and files
        self.ligand_dir = os.path.join(self.dataset_dir, 'ligand_files')
        self.model_dir  = os.path.join(self.dataset_dir, 'modelled_structures')

        if not os.path.exists(self.model_dir): os.mkdir(self.model_dir)

        # The most recent model of the protein in the pandda maps
        self.fitted_link = os.path.join(self.model_dir, 'fitted-current.pdb')
        # Unfitted model of the protein
        self.unfitted_model = os.path.join(self.dataset_dir, '{!s}-aligned.pdb'.format(dtag))
        # Maps
        self.observed_map   = os.path.join(self.dataset_dir, '{!s}-observed.ccp4'.format(dtag))
        self.z_map          = os.path.join(self.dataset_dir, '{!s}-z_map_adjusted_normalised.ccp4'.format(dtag))
        self.mean_diff_map  = os.path.join(self.dataset_dir, '{!s}-mean_diff.ccp4'.format(dtag))
        try: self.occupancy_map = glob.glob(os.path.join(self.dataset_dir, '{!s}-event_{!s}_occupancy_*_map.ccp4'.format(dtag, self.event_idx)))[0]
        except: raise

        # Ligand Files
        lig_files = glob.glob(os.path.join(self.ligand_dir, '*'))
        self.lig_pdbs = [f for f in lig_files if f.endswith('.pdb')]
        self.lig_cifs = [f for f in lig_files if f.endswith('.cif')]

    def find_current_fitted_model(self):
        """Get the most recent saved model of this protein"""
        fitted_outputs = sorted(glob.glob(os.path.join(self.model_dir, 'fitted-v*')))
        if fitted_outputs: return fitted_outputs[-1]
        else:              return None

    def find_new_fitted_model(self):
        """Get the path for the next model"""
        current = self.find_current_fitted_model()
        print 'CURRENT: {!s}'.format(current)
        if current: last_idx = int(current.replace('.pdb','')[-4:])
        else:       last_idx = 0
        new_fitted = 'fitted-v{:04d}.pdb'.format(last_idx+1)
        return os.path.join(self.model_dir, new_fitted)

    def write_fitted_model(self, protein_obj):
        new_fitted = self.find_new_fitted_model()
        # Save the fitted model to this file
        write_pdb_file(protein_obj, new_fitted)
        print '======================>'
        print 'SAVED PROTEIN MOL {!s} TO {!s}'.format(protein_obj, new_fitted)
        print '======================>'
        if os.path.exists(self.fitted_link):
            assert os.path.islink(self.fitted_link), 'FILE IS NOT A LINK'
            # Delete the symbolic link
            os.unlink(self.fitted_link)
        print 'Linking {!s} -> {!s}'.format(new_fitted, self.fitted_link)
        # Create new link the most recent file
        os.symlink(os.path.basename(new_fitted), self.fitted_link)

#=========================================================================

class PanddaSiteTracker(object):
    events = None
    # Sites count from 1
    site_val = 1
    site_tot = 1
    # Rank counts from 1
    rank_val = 1
    rank_idx = 0
    rank_tot = 1

    def __init__(self, csv, top_dir):
        self.top_dir = top_dir
        self.parse_csv(csv)

    def parse_csv(self, csv):

        self.events = pandas.read_csv(csv, sep=',', dtype={'dtag':str})
        self.events = self.events.set_index(['dtag','event_idx'])

        print '====================>>>'
        print self.events
        print '====================>>>'
        current = self.get_new_current_event()
        print current
        print '====================>>>'

        self.update()

    def update(self):
        self.site_tot = len(set(self.events['site_idx']))
        self.rank_tot = len(self.events)

    def get_new_current_event(self):
        """Get an `Event` object for the current event"""
        curr_event = self.events.iloc[self.rank_idx]
        self.site_val = int(curr_event['site_idx'])
        print '\n\nCurrent Event:\n\n{!s}\n\n'.format(curr_event)
        return PanddaEvent(rank=self.rank_val, info=curr_event, top_dir=self.top_dir)

    #-------------------------------------------------------------------------

    def get_next(self):
        if self.rank_val == self.rank_tot: return None
        self.rank_idx += 1
        self.rank_val = self.rank_idx + 1
        return self.get_new_current_event()

    def get_prev(self):
        if self.rank_val == 1: return None
        self.rank_idx -= 1
        self.rank_val = self.rank_idx + 1
        return self.get_new_current_event()

    def get_prev_site(self):
        if sel.site_val == 1:  return None
        self.rank_idx = self.events.index.get_loc(self.events[self.events['site_idx']==self.site_val-1].index[0])
        self.rank_val = self.rank_idx + 1
        return self.get_new_current_event()

    def get_next_site(self):
        if curr_site == self.site_tot: return None
        self.rank_idx = self.events.index.get_loc(self.events[self.events['site_idx']==self.site_val+1].index[0])
        self.rank_val = self.rank_idx + 1
        return self.get_new_current_event()

#=========================================================================

class PanddaInspector(object):

    """Main Object in pandda.inspect"""

    def __init__(self, csv, top_dir):

        # List of events from pandda.analyse
        self.site_list = PanddaSiteTracker(csv=csv, top_dir=top_dir)
        # Handling of coot commands
        self.coot = PanddaMolHandler(parent=self)

        # Working Directory
        self.top_dir = top_dir
        # Output csv from pandda.inspect
        self.output_csv = os.path.join(self.top_dir, 'pandda_inspect.csv')

        # Load previous data or create new table from site list, so we can record inspection data
        self.initialise_output_table(in_csv=csv, out_csv=self.output_csv)

        # Object to hold the currently loaded event
        self.current_event = None

    def start_gui(self):
        self.gui = PanddaGUI(parent=self)
        self.gui.launch()

    def update_gui(self):

        # Update global position
        self.gui.labels['site_val'].set_label(str(self.site_list.site_val))
        self.gui.labels['site_tot'].set_label(str(self.site_list.site_tot))
        self.gui.labels['rank_val'].set_label(str(self.site_list.rank_val))
        self.gui.labels['rank_tot'].set_label(str(self.site_list.rank_tot))
        # Update current event information
        self.gui.labels['dtag'].set_label(str(self.current_event.dtag))
        self.gui.labels['idx'].set_label(str(self.current_event.event_idx))
        self.gui.labels['zpeak'].set_label(str(self.current_event.z_peak))
        self.gui.labels['csize'].set_label(str(self.current_event.cluster_size))

    #-------------------------------------------------------------------------

    def load_next_event(self, skip_unmodelled=False):
        # Get the next event from the site list
        new_event = self.site_list.get_next()
        # Loop until we get a modelled one (if requested)
        if skip_unmodelled and (not os.path.exists(new_event.fitted_link)):
            self.load_next_event(); return
        # Load and Save the event
        self.current_event = self.coot.load_event(e=new_event)
        # Update the gui
        self.update_gui()

    def load_prev_event(self, skip_unmodelled=False):
        # Get the next event from the site list
        new_event = self.site_list.get_prev()
        # Loop until we get a modelled one (if requested)
        if skip_unmodelled and (not os.path.exists(new_event.fitted_link)):
            self.load_prev_event(); return
        # Load and Save the event
        self.current_event = self.coot.load_event(e=new_event)
        # Update the gui
        self.update_gui()

    def refresh_event(self):
        # Get new event instance
        new_event = self.site_list.get_new_current_event()
        # Load and store the event
        self.current_event = self.coot.load_event(e=new_event)
        # Update the gui
        self.update_gui()

    #-------------------------------------------------------------------------

    def save_and_load_next(self, skip_unmodelled=False):
        self.save_current()
        self.load_next_event(skip_unmodelled=skip_unmodelled)

    def save_and_load_prev(self, skip_unmodelled=False):
        self.save_current()
        self.load_prev_event(skip_unmodelled=skip_unmodelled)

    def save_current(self):
        self.current_event.write_fitted_model(protein_obj=self.coot.open_mols['p'])
        self.write_output_csv()

    #-------------------------------------------------------------------------

    def set_log_value(self, col, value):
        self.log_table.set_value(index=self.current_event.index, col=col, value=value)
        print '====================>>>'
        print 'LOG TABLE UPDATED:'
        print '====================>>>'
        print self.log_table.loc[self.current_event.index]
        print '====================>>>'
        self.write_output_csv()

    #-------------------------------------------------------------------------

    def initialise_output_table(self, in_csv, out_csv):

        if os.path.exists(out_csv):
            # Output csv already exists from previous run -- reload
            self.log_table = pandas.read_csv(out_csv, sep=',', dtype={'dtag':str})
        else:
            # Create new table from input csv
            self.log_table = pandas.read_csv(in_csv, sep=',', dtype={'dtag':str})
            # Create new columns
            self.log_table['Interesting'] = False
            self.log_table['Ligand Placed'] = False

        # Set the index
        self.log_table = self.log_table.set_index(['dtag','event_idx'])

        # Save Table
        self.write_output_csv()

    def write_output_csv(self):
        self.log_table.to_csv(self.output_csv)

#=========================================================================

class PanddaMolHandler(object):

    """Handles loaded Pandda Models (contains most of the coot functions)"""

    def __init__(self, parent):
        self.parent = parent
        self.open_mols = {}

    def close_all(self):
        for mol in molecule_number_list():
            close_molecule(mol)
        self.open_mols = {}

    def merge_ligand_with_protein(self):
        try: merge_molecules([self.open_mols['l']], self.open_mols['p'])
        except Exception as err: print err

    def move_ligand_here(self):
        try: move_molecule_to_screen_centre(self.open_mols['l'])
        except Exception as err: print err

    def load_event(self, e, close_all=True):
        """Load up all of the maps for an event"""

        # Close 'n' Load
        if close_all: self.close_all()

        # Re-centre camera
        set_rotation_centre(*e.ref_coords)

        # Load fitted version if it exists
        if os.path.exists(e.fitted_link): p = read_pdb(e.fitted_link)
        else:                             p = read_pdb(e.unfitted_model)

        # Load fitted version if it exists
        s = handle_read_ccp4_map(e.observed_map, 0)
        set_last_map_contour_level(1)
        set_map_displayed(s, 0)

        # Load fitted version if it exists
        z = handle_read_ccp4_map(e.z_map, 1)
        set_last_map_contour_level(3)
        set_map_displayed(z, 1)

        # Mean-Difference Map
        d = handle_read_ccp4_map(e.mean_diff_map, 1)
        set_last_map_contour_level_by_sigma(3)
        set_map_displayed(d, 0)

        # Occupancy Map
#        try:
        o = handle_read_ccp4_map(e.occupancy_map, 0)
        set_last_map_contour_level(1)
        set_map_displayed(o, 1)
#        except: o = z

        # Symmetry contacts
        r = read_pdb(e.reference_symmetry)
        set_mol_displayed(r, 0)

        # Ligand Files
        if (len(e.lig_pdbs) == 1) and (len(e.lig_cifs) == 1):
            l_dict = read_cif_dictionary(e.lig_cifs[0])
            l = handle_read_draw_molecule_and_move_molecule_here(e.lig_pdbs[0])
            if os.path.exists(e.fitted_link): set_mol_displayed(l, 0)
            else:                             set_mol_displayed(l, 1)

        # More Settings
        set_scrollable_map(o)
        set_imol_refinement_map(o)

        # Save mol numbers
        self.open_mols['p'] = p
        self.open_mols['s'] = s
        self.open_mols['z'] = z
        self.open_mols['d'] = d
        self.open_mols['o'] = o
        self.open_mols['r'] = r
        self.open_mols['l'] = l

        return e

class GenericGUI(object):

    """Generic GUI Class"""

    def generic_label_value_pair(self, label, value):
        gtk_label = gtk.Label(label)
        gtk_value = gtk.Label(value)
        gtk_box = gtk.EventBox()
        gtk_box.add(gtk_value)
        hbox = gtk.HBox()
        hbox.add(gtk_label)
        hbox.add(gtk_box)
        frame = gtk.Frame()
        frame.add(hbox)
        return (gtk_value, frame)

class PanddaGUI(GenericGUI):

    """GUI Class for pandda.inspect"""

    def __init__(self, parent):

        self.parent = parent

        # GUI objects
        self.labels = {}
        self.buttons = None
        self.objects = {}

    def launch(self):
        """Launch GUI window"""

        # Create main window
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.connect("delete_event", gtk.main_quit)
        self.window.set_border_width(10)
        self.window.set_default_size(600, 400)
        self.window.set_title("PANDDA inspect")

        # Main VBox object for the window
        main_vbox = gtk.VBox()
        main_vbox.set_spacing(5)

        # Create summary table at the top of the window
        self.progress_table = self.new_progress_table()
        self.event_table    = self.new_event_table()
        # Create own vbox
        label_vbox = gtk.VBox()
        label_vbox.set_spacing(5)
        frame = gtk.Frame(); frame.add(self.progress_table)
        label_vbox.add(frame)
        frame = gtk.Frame(); frame.add(self.event_table)
        label_vbox.add(frame)
        # Add vbox to main vbox
        main_vbox.add(label_vbox)

        # Create buttons at the bottom of the window (button_vbox is populated)
        button_vbox, self.buttons = self.main_buttons()
        # Add vbox to main vbox
        main_vbox.add(button_vbox)

        # Link the buttons to the Inspector
        self.link_buttons()

        # Finally, show the window
        self.window.add(main_vbox)
        self.window.show_all()

        return self

    def main_buttons(self):

        button_vbox = gtk.VBox()
        button_vbox.set_spacing(5)
        button_dict = {}

        # Navigation Buttons
        nav_hbox = gtk.HBox()
        # ---
        prev = gtk.Button(label="<<< Prev <<<")
        button_dict['prev'] = prev
        nav_hbox.add(prev)
        # ---
        skip = gtk.Button(label=">>> Next (Skip) >>>")
        button_dict['skip'] = skip
        nav_hbox.add(skip)
        # ---
        next = gtk.Button(label=">>> Next (+Save) >>>")
        button_dict['next'] = next
        nav_hbox.add(next)
        # ---
        button_vbox.add(nav_hbox)

        # Data Buttons
        data_hbox = gtk.HBox()
        vbox_1 = gtk.VBox()
        vbox_2 = gtk.VBox()
        # ---
        b = gtk.Button(label="Mark Interesting")
        button_dict['tp'] = b
        vbox_1.add(b)
        # ---
        b = gtk.Button(label="Mark Not Interesting")
        button_dict['fp'] = b
        vbox_1.add(b)
        # ---
        b = gtk.Button(label="Ligand Placed")
        button_dict['placed'] = b
        vbox_2.add(b)
        # ---
        data_hbox.add(vbox_1)
        data_hbox.add(vbox_2)
        button_vbox.add(data_hbox)

        # Helper Buttons
        helper_hbox = gtk.HBox()
        # ---
        merge = gtk.Button(label="Merge Ligand With Model")
        button_dict['merge'] = merge
        helper_hbox.add(merge)
        # ---
        move = gtk.Button(label="Move New Ligand Here")
        button_dict['move'] = move
        helper_hbox.add(move)
        # ---
        save = gtk.Button(label="Save Model")
        button_dict['save'] = save
        helper_hbox.add(save)
        # ---
        button_vbox.add(helper_hbox)

        return button_vbox, button_dict

    def link_buttons(self):
        """Link the buttons in the GUI to functions in PanddaInspector"""

        # Navigation buttons
        self.buttons['next'].connect("clicked", lambda x: self.parent.save_and_load_next())
        self.buttons['prev'].connect("clicked", lambda x: self.parent.save_and_load_prev())
        self.buttons['skip'].connect("clicked", lambda x: self.parent.load_next_event())

        # Structure Buttons
        self.buttons['save'].connect("clicked",  lambda x: self.parent.save_current())
        self.buttons['merge'].connect("clicked", lambda x: self.parent.coot.merge_ligand_with_protein())
        self.buttons['move'].connect("clicked",  lambda x: self.parent.coot.move_ligand_here())

        # Meta Recording buttons
        self.buttons['placed'].connect("clicked", lambda x: self.parent.set_log_value(col='Ligand Placed', value='True'))
        self.buttons['tp'].connect("clicked", lambda x: self.parent.set_log_value(col='Interesting', value=True))
        self.buttons['fp'].connect("clicked", lambda x: self.parent.set_log_value(col='Interesting', value=False))

    def new_event_table(self):

        self.labels['dtag'],  tag_pair   = self.generic_label_value_pair( label='Dataset', value='None' )
        self.labels['idx'],   idx_pair   = self.generic_label_value_pair( label='Event #', value=1 )
        self.labels['zpeak'], zpeak_pair = self.generic_label_value_pair( label='Z-Peak',  value=0 )
        self.labels['csize'], csize_pair = self.generic_label_value_pair( label='Z-Size',  value=0 )

        table = gtk.Table(rows=3, columns=4)
        title = gtk.Label('Event Information:')
        title.set_justify(gtk.JUSTIFY_LEFT)
        table.attach(title, 0,4,0,1)
        table.attach(tag_pair,   0,2,1,2)
        table.attach(idx_pair,   2,4,1,2)
        table.attach(zpeak_pair, 0,2,2,3)
        table.attach(csize_pair, 2,4,2,3)

        return table

    def new_progress_table(self):

        self.labels['rank_val'], rank_val_pair = self.generic_label_value_pair( label='Event', value=0 )
        self.labels['rank_tot'], rank_tot_pair = self.generic_label_value_pair( label='Total',    value=0 )
        self.labels['site_val'], site_val_pair = self.generic_label_value_pair( label='Site',  value=0 )
        self.labels['site_tot'], site_tot_pair = self.generic_label_value_pair( label='Total',    value=0 )

        table = gtk.Table(rows=2, columns=8)
        title = gtk.Label('Inspection Event/Site Progress:')
        title.set_justify(gtk.JUSTIFY_LEFT)
        table.attach(title, 0,4,0,2)
        table.attach(rank_val_pair, 4,6,0,1)
        table.attach(rank_tot_pair, 6,8,0,1)
        table.attach(site_val_pair, 4,6,1,2)
        table.attach(site_tot_pair, 6,8,1,2)

        return table

#############################################################################################

#try:
#    add_key_binding("next blob", "s", lambda: inspector.save_and_load_next())
#    add_key_binding("prev blob", "a", lambda: inspector.save_and_load_prev())
#    add_key_binding("skip blob", "v", lambda: inspector.load_next_event())
#
#except Exception as err:
#    print err

#############################################################################################



#############################################################################################
#
# INITIALISE THE MAIN INSPECTOR OBJECT
#
#############################################################################################

work_dir = os.getcwd()
hit_list = os.path.join(work_dir, 'analyses', 'event_info.csv')

inspector = PanddaInspector(csv=hit_list, top_dir=work_dir)
inspector.start_gui()
inspector.refresh_event()


