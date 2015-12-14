
import os, sys, glob, copy, itertools
import gtk
import pandas

from PANDDAs import graphs, inspect_html
from PANDDAs.constants import PanddaAnalyserFilenames, PanddaDatasetFilenames, PanddaDatasetPNGFilenames

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
        # Event Info
        self.est_occ = round(info['est_occupancy'], 2)
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
        # Find the directory for the dataset
        self.dataset_dir  = glob.glob(os.path.join(top_dir, 'interesting_datasets/{!s}'.format(dtag)))
        assert len(self.dataset_dir) == 1, 'Non-unique dataset directory found: {!s}'.format(dtag)
        self.dataset_dir  = self.dataset_dir[0]

        # Identify dataset subdirectories and files
        self.ligand_dir = os.path.join(self.dataset_dir, 'ligand_files')
        self.model_dir  = os.path.join(self.dataset_dir, 'modelled_structures')

        if not os.path.exists(self.model_dir): os.mkdir(self.model_dir)

        # The most recent model of the protein in the pandda maps
        self.fitted_link        = os.path.join(self.model_dir,      PanddaDatasetFilenames.modelled_structure.format(dtag)  )
        # Unfitted model of the protein
        self.unfitted_model     = os.path.join(self.dataset_dir,    PanddaDatasetFilenames.aligned_structure.format(dtag)   )
        # Maps
        self.observed_map       = os.path.join(self.dataset_dir,    PanddaDatasetFilenames.sampled_map.format(dtag)         )
        self.z_map              = os.path.join(self.dataset_dir,    PanddaDatasetFilenames.z_map.format(dtag)               )
        self.mean_diff_map      = os.path.join(self.dataset_dir,    PanddaDatasetFilenames.mean_diff_map.format(dtag)       )
        self.occupancy_map      = os.path.join(self.dataset_dir,    PanddaDatasetFilenames.occupancy_map.format(dtag, self.event_idx, self.est_occ))
        self.dataset_symmetry   = os.path.join(self.dataset_dir,    PanddaDatasetFilenames.symmetry_copies.format(dtag)     )

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
    # Sites
    site_val = 1 # Site Value           1 -> m
    site_tot = 1 # Number of Sites      m
    # Ranks (Events)
    rank_idx = 0 # Indexing position    0 -> n-1
    rank_val = 1 # Rank Value           1 -> n
    rank_tot = 1 # Number of Events     n

    def __init__(self, parent, csv, top_dir):
        self.parent  = parent
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
        """Update position of the site list"""
        self.rank_val = self.rank_idx + 1
        self.rank_tot = len(self.events)
        self.site_val = int(self.events.iloc[self.rank_idx]['site_idx'])
        self.site_tot = len(set(self.events['site_idx']))

    def get_new_current_event(self):
        """Get an `Event` object for the current event"""
        self.update() # Ensure that we're up-to-date
        curr_event = self.events.iloc[self.rank_idx]
#        print '\n\nCurrent Event:\n\n{!s}\n\n'.format(curr_event)
        return PanddaEvent(rank=self.rank_val, info=curr_event, top_dir=self.top_dir)

    #-------------------------------------------------------------------------

    def at_first_event(self):
        self.update()
        return self.rank_val == 1

    def at_last_event(self):
        self.update()
        return self.rank_val == self.rank_tot

    #-------------------------------------------------------------------------

    def get_next(self):
        if self.rank_idx == self.rank_tot - 1:
            self.rank_idx = 0
#            return None
        else:
            self.rank_idx += 1
        return self.get_new_current_event()

    def get_prev(self):
        if self.rank_idx == 0:
            self.rank_idx = self.rank_tot - 1
#            return None
        else:
            self.rank_idx -= 1
        return self.get_new_current_event()

    def get_next_site(self):
        if self.site_val == self.site_tot:
            new_site_val = 1
#            return None
        else:
            new_site_val = self.site_val + 1
        self.rank_idx = self.events.index.get_loc(self.events[self.events['site_idx']==new_site_val].index[0])
        return self.get_new_current_event()

    def get_prev_site(self):
        if self.site_val == 1:
            new_site_val = self.site_tot
#            return None
        else:
            new_site_val = self.site_val - 1
        self.rank_idx = self.events.index.get_loc(self.events[self.events['site_idx']==new_site_val].index[0])
        return self.get_new_current_event()

#=========================================================================

class PanddaInspector(object):

    """Main Object in pandda.inspect"""

    def __init__(self, csv, top_dir):

        # List of events from pandda.analyse
        self.site_list = PanddaSiteTracker(parent=self, csv=csv, top_dir=top_dir)
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

    def update_html(self):
        """Update various parts of the class - including output graphs"""

        # Plot output graph of site list
        plot_vals = self.log_table['z_peak']
        view_vals = self.log_table['Viewed']
        modl_vals = self.log_table['Ligand Placed']
        colr_vals = ['limegreen' if m else 'red' if v else 'blue' for m,v in zip(modl_vals,view_vals)]
        site_idxs = self.log_table['site_idx']
        groups = [list(g[1]) for g in itertools.groupby(range(len(site_idxs)), key=lambda i: site_idxs[i])]
        graphs.multiple_bar_plot(   f_name      = os.path.join(self.top_dir, 'results_summaries', 'pandda_inspect_sites.png'),
                                    plot_vals   = [[plot_vals[i] for i in g] for g in groups],
#                                    colour_bool = [[view_vals[i] for i in g] for g in groups]   )
                                    colour_vals = [[colr_vals[i] for i in g] for g in groups]   )
        # Write output html
        inspect_html.write_inspect_html(out_dir=self.top_dir, inspector=self)

    def update_gui(self):

        # Update global position
        self.gui.labels['site_val'].set_label(str(self.site_list.site_val))
        self.gui.labels['site_tot'].set_label(str(self.site_list.site_tot))
        self.gui.labels['rank_val'].set_label(str(self.site_list.rank_val))
        self.gui.labels['rank_tot'].set_label(str(self.site_list.rank_tot))
        # Update current event information
        self.gui.labels['dtag'].set_label(str(self.current_event.dtag))
        self.gui.labels['e_idx'].set_label(str(self.current_event.event_idx))
        self.gui.labels['e_occ'].set_label(str(self.current_event.est_occ))
        self.gui.labels['zpeak'].set_label(str(round(self.current_event.z_peak,3)))
        self.gui.labels['csize'].set_label(str(self.current_event.cluster_size))
        # Reset the comment box
        self.gui.objects['comment text'].set_text(str(self.get_log_value('Comment')))

    def raise_gui_error_and_reset(self, error_msg):
        self.gui.error_msg(msg=error_msg)

    #-------------------------------------------------------------------------

    def save_current(self):
        self.current_event.write_fitted_model(protein_obj=self.coot.open_mols['p'])
        self.write_output_csv()

    def load_new_event(self, new_event):
        if new_event is None:
            # TODO DO SOMETHING BETTER HERE TODO
            # Instead of loading new event, inform using the same model
            self.raise_gui_error_and_reset(error_msg='Reloading Same MODEL (Unchanged)')
        else:
            self.current_event = self.coot.load_event(e=new_event)
        self.update_gui()
        self.print_current_log_values()

    def refresh_event(self):
        self.load_new_event(new_event=self.site_list.get_new_current_event())

    def load_next_event(self, skip_unmodelled=False, skip_viewed=False):
        if not self.check_valid_button_click(skip_unmodelled, skip_viewed): return
        new_event = self.site_list.get_next()
        # Loop until we get a modelled one (if requested)
        if (new_event is not None) and skip_unmodelled and (not os.path.exists(new_event.fitted_link)):
#            if self.site_list.at_first_event():
            if self.current_event.index == new_event.index:
                self.raise_gui_error_and_reset(error_msg='No more unmodelled datasets')
            else:
                print 'SKIPPING UNMODELLED: {} - {}'.format(new_event.dtag, new_event.event_idx)
                self.load_next_event(skip_unmodelled=skip_unmodelled, skip_viewed=skip_viewed)
            return
        # Loop until we get an unviewed one (if requested)
        if (new_event is not None) and skip_viewed and (self.log_table.get_value(index=new_event.index, col='Viewed') == True):
            print 'SKIPPING VIEWED: {} - {}'.format(new_event.dtag, new_event.event_idx)
            self.load_next_event(skip_unmodelled=skip_unmodelled, skip_viewed=skip_viewed)
            return
        # Actually load the event
        self.load_new_event(new_event=new_event)

    def load_prev_event(self, skip_unmodelled=False, skip_viewed=False):
        if not self.check_valid_button_click(skip_unmodelled, skip_viewed): return
        new_event = self.site_list.get_prev()
#        # Loop until we get a modelled one (if requested)
#        if (new_event is not None) and skip_unmodelled and (not os.path.exists(new_event.fitted_link)):
#            print 'SKIPPING UNMODELLED: {} - {}'.format(new_event.dtag, new_event.event_idx)
#            self.load_prev_event(skip_unmodelled=skip_unmodelled, skip_viewed=skip_viewed)
#            return
#        # Loop until we get an unviewed one (if requested)
#        if (new_event is not None) and skip_viewed and (self.log_table.get_value(index=new_event.index, col='Viewed') == True):
#            print 'SKIPPING VIEWED: {} - {}'.format(new_event.dtag, new_event.event_idx)
#            self.load_prev_event(skip_unmodelled=skip_unmodelled, skip_viewed=skip_viewed)
#            return
        # Actually load the event
        self.load_new_event(new_event=new_event)

    def load_next_site(self):
        self.load_new_event(new_event=self.site_list.get_next_site())

    def load_prev_site(self):
        self.load_new_event(new_event=self.site_list.get_prev_site())

    #-------------------------------------------------------------------------

    def check_valid_button_click(self, skip_unmodelled=None, skip_viewed=None):
        if skip_unmodelled == True:
            pass
        if skip_viewed == True:
            if sum(self.log_table['Viewed']) == len(self.log_table['Viewed']):
                self.raise_gui_error_and_reset(error_msg='All models have been viewed')
                return False
        return True
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
            self.log_table['Ligand Confidence'] = 'None'
            self.log_table['Comment'] = 'None'
            self.log_table['Viewed'] = False
        # Set the index
        self.log_table = self.log_table.set_index(['dtag','event_idx'])
        # Save Table
        self.write_output_csv()

    def get_log_value(self, col):
        return self.log_table.get_value(index=self.current_event.index, col=col)

    def set_log_value(self, col, value):
        self.log_table.set_value(index=self.current_event.index, col=col, value=value)
        print '====================>>>'
        print 'LOG TABLE UPDATED:'
        self.print_current_log_values()
        self.write_output_csv()

    def print_current_log_values(self):
        print '====================>>>'
        print 'Current Event:'
        print self.log_table.loc[self.current_event.index]
        print '====================>>>'

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
        o = handle_read_ccp4_map(e.occupancy_map, 0)
        set_last_map_contour_level(1)
        set_map_displayed(o, 1)

        # Symmetry contacts
        r = read_pdb(e.dataset_symmetry)
        set_mol_displayed(r, 0)

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

        # Ligand Files
        if (len(e.lig_pdbs) == 1) and (len(e.lig_cifs) == 1):
            l_dict = read_cif_dictionary(e.lig_cifs[0])
            l = handle_read_draw_molecule_and_move_molecule_here(e.lig_pdbs[0])
            if os.path.exists(e.fitted_link): set_mol_displayed(l, 0)
            else:                             set_mol_displayed(l, 1)
            self.open_mols['l'] = l

        return e

class PanddaGUI(object):

    """GUI Class for pandda.inspect"""

    def __init__(self, parent):

        self.parent = parent

        # GUI objects
        self.labels = {}
        self.buttons = {}
        self.objects = {}

    def on_destroy(self,  widget=None, *data):
        gtk.main_quit()

    def error_msg(self, msg):
        """Display an error window"""
        d = gtk.MessageDialog(  parent  = self.window,
                                flags   = gtk.DIALOG_DESTROY_WITH_PARENT,
                                type    = gtk.MESSAGE_INFO,
                                buttons = gtk.BUTTONS_CLOSE,
                                message_format = msg )
        d.run()
        d.destroy()

    def launch(self):
        """Launch GUI window"""

        # Create main window
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.connect("delete_event", self.on_destroy)
        self.window.connect("destroy_event", self.on_destroy)
        self.window.set_border_width(10)
        self.window.set_default_size(600, 400)
        self.window.set_title("PANDDA inspect")

        # Main VBox object for the window
        main_vbox = gtk.VBox(spacing=5)
        # -----------------------------------------------------
        hbox = gtk.HBox(spacing=5)
        # Create progress summary table at the top of the window
        self.progress_table = self._progress_table()
        frame = gtk.Frame(); frame.add(self.progress_table)
        hbox.pack_start(frame)
        # Create buttons to navigate between datasets
        nav_buttons = self._navi_buttons_2()
        hbox.pack_start(nav_buttons)
        # Add to main vbox
        main_vbox.pack_start(hbox)
        # -----------------------------------------------------
        # Create buttons to navigate between datasets
        nav_buttons = self._navi_buttons_1()
        frame = gtk.Frame(); frame.add(nav_buttons)
        main_vbox.pack_start(frame)
        # -----------------------------------------------------
        hbox = gtk.HBox(homogeneous=False, spacing=5)
        # Create event summary table at the top of the window
        self.event_info_table = self._event_info_table()
        frame = gtk.Frame(); frame.add(self.event_info_table)
        hbox.pack_start(frame)
        # Create buttons to control the ligand
        lig_buttons = self._ligand_buttons()
        frame = gtk.Frame(); frame.add(lig_buttons)
        hbox.pack_start(frame)
        # Add to main vbox
        main_vbox.pack_start(hbox)
        # -----------------------------------------------------
        # Create buttones to record meta about the event
        rec_buttons = self._record_buttons()
        frame = gtk.Frame(); frame.add(rec_buttons)
        main_vbox.pack_start(frame)

        # Link the buttons to the Inspector
        self.link_buttons()

        # Finally, show the window
        self.window.add(main_vbox)
        self.window.show_all()

        return self

    def link_buttons(self):
        """Link the buttons in the GUI to functions in PanddaInspector"""

        # Navigation buttons
        self.buttons['next'].connect("clicked", lambda x: [self.store(), self.parent.save_current(), self.parent.load_next_event()])
        self.buttons['prev'].connect("clicked", lambda x: [self.store(), self.parent.load_prev_event()])
        self.buttons['skip'].connect("clicked", lambda x: [self.store(), self.parent.load_next_event()])

        self.buttons['next-unviewed'].connect("clicked", lambda x: [self.store(), self.parent.load_next_event(skip_viewed=True)])
        self.buttons['next-modelled'].connect("clicked", lambda x: [self.store(), self.parent.load_next_event(skip_unmodelled=True)])
        self.buttons['next-site'].connect("clicked", lambda x: [self.store(), self.parent.load_next_site()])
        self.buttons['prev-site'].connect("clicked", lambda x: [self.store(), self.parent.load_prev_site()])

        # Structure Buttons
        self.buttons['save'].connect("clicked",  lambda x: self.parent.save_current())
        self.buttons['merge'].connect("clicked", lambda x: self.parent.coot.merge_ligand_with_protein())
        self.buttons['move'].connect("clicked",  lambda x: self.parent.coot.move_ligand_here())

        # Meta Recording buttons
        self.buttons['tp'].connect("clicked",         lambda x: self.parent.set_log_value(col='Interesting', value=True))
        self.buttons['fp'].connect("clicked",         lambda x: self.parent.set_log_value(col='Interesting', value=False))
        self.buttons['high conf'].connect("clicked",  lambda x: self.parent.set_log_value(col='Ligand Confidence', value='High'))
        self.buttons['med conf'].connect("clicked",   lambda x: self.parent.set_log_value(col='Ligand Confidence', value='Medium'))
        self.buttons['low conf'].connect("clicked",   lambda x: self.parent.set_log_value(col='Ligand Confidence', value='Low'))
        self.buttons['placed'].connect("clicked",     lambda x: self.parent.set_log_value(col='Ligand Placed', value=True))
        self.buttons['not placed'].connect("clicked", lambda x: self.parent.set_log_value(col='Ligand Placed', value=False))

    def store(self):
        """Record information from the gui to the pandas table in the main object"""
        self.parent.set_log_value(col='Comment', value=self.objects['comment text'].get_text())
        self.parent.set_log_value(col='Viewed', value=True)
        self.parent.update_html()

    def _navi_buttons_1(self):
        box = gtk.HBox(homogeneous=False, spacing=2)
        box.set_border_width(3)
        # ---
        b = gtk.Button(label="<<< Prev (Skip) <<<")
        self.buttons['prev'] = b
        box.pack_start(b)
        # ---
        b = gtk.Button(label=">>> Next (Skip) >>>")
        self.buttons['skip'] = b
        box.pack_start(b)
        # ---
        b = gtk.VSeparator()
        box.pack_start(b, expand=False, padding=5)
        # ---
        b = gtk.Button(label=">>> Next (Save) >>>")
        self.buttons['next'] = b
        box.pack_start(b)
        # ---
        return box

    def _navi_buttons_2(self):
        main_box  = gtk.VBox(homogeneous=True, spacing=5)
        # ---
        box1 = gtk.HBox(homogeneous=True, spacing=2)
        box1.set_border_width(3)
        frame = gtk.Frame(); frame.add(box1)
        main_box.pack_start(frame)
        # ---
        box2 = gtk.HBox(homogeneous=True, spacing=2)
        box2.set_border_width(3)
        frame = gtk.Frame(); frame.add(box2)
        main_box.pack_start(frame)
        # ---
        b = gtk.Button(label="<<< Prev Site <<<")
        self.buttons['prev-site'] = b
        box1.add(b)
        # ---
        b = gtk.Button(label=">>> Next Site >>>")
        self.buttons['next-site'] = b
        box1.add(b)
        # ---
        b = gtk.Button(label=">>> Next Unviewed >>>")
        self.buttons['next-unviewed'] = b
        box2.add(b)
        # ---
        b = gtk.Button(label=">>> Next Modelled >>>")
        self.buttons['next-modelled'] = b
        box2.add(b)
        # ---
        return main_box

    def _ligand_buttons(self):
        box = gtk.HBox(homogeneous=True, spacing=2)
        box.set_border_width(3)
        # ---
        b = gtk.Button(label="Merge Ligand With Model")
        b.child.set_line_wrap(True)
        b.child.props.width_chars = 10
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        self.buttons['merge'] = b
        box.add(b)
        # ---
        b = gtk.Button(label="Move New Ligand Here")
        b.child.set_line_wrap(True)
        b.child.props.width_chars = 10
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        self.buttons['move'] = b
        box.add(b)
        # ---
        b = gtk.Button(label="Save Model")
        self.buttons['save'] = b
        box.add(b)
        # ---
        return box

    def _record_buttons(self):
        # ---------------------------------------------
        hbox_1 = gtk.HBox(homogeneous=True, spacing=5)
        hbox_1.set_border_width(3)
        vbox_1_1 = gtk.VBox(homogeneous=True, spacing=2)
        vbox_1_2 = gtk.VBox(homogeneous=True, spacing=2)
        vbox_1_3 = gtk.VBox(homogeneous=True, spacing=2)
        hbox_1.add(vbox_1_1); hbox_1.add(vbox_1_2); hbox_1.add(vbox_1_3)
        # ---
        b = gtk.Button(label="Mark Interesting")
        self.buttons['tp'] = b
        vbox_1_1.add(b)
        # ---
        b = gtk.Button(label="Mark Not Interesting")
        self.buttons['fp'] = b
        vbox_1_1.add(b)
        # ---
        b = gtk.Button(label="High Confidence")
        self.buttons['high conf'] = b
        vbox_1_2.add(b)
        # ---
        b = gtk.Button(label="Medium Confidence")
        self.buttons['med conf'] = b
        vbox_1_2.add(b)
        # ---
        b = gtk.Button(label="Low Confidence")
        self.buttons['low conf'] = b
        vbox_1_2.add(b)
        # ---
        b = gtk.Button(label="Ligand Placed")
        self.buttons['placed'] = b
        vbox_1_3.add(b)
        # ---
        b = gtk.Button(label="No Ligand Placed")
        self.buttons['not placed'] = b
        vbox_1_3.add(b)
        # ---------------------------------------------
        hbox_2 = gtk.HBox(homogeneous=False, spacing=5)
        # ---
        l = gtk.Label('Comment:')
        hbox_2.add(l)
        # ---
        e = gtk.Entry(max=200)
        self.objects['comment text'] = e
        hbox_2.add(e)
        # ---------------------------------------------
        vbox_main = gtk.VBox(spacing=5)
        vbox_main.pack_start(hbox_1)
        vbox_main.pack_start(hbox_2)

        return vbox_main

    def _event_info_table(self):

        # First Column
        vbox_1 = gtk.VBox()
        vbox_1.set_border_width(3)

        # Create title
        title = gtk.Label('Event Information:')
        title.set_justify(gtk.JUSTIFY_LEFT)
        frame = gtk.Frame(); frame.add(title)
        # Add to first column
        vbox_1.pack_start(frame)

        # Dataset Name
        gtk_label = gtk.Label('Dataset')
        gtk_value = gtk.Label('None')
        gtk_box = gtk.EventBox(); gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True); hbox.add(gtk_label); hbox.add(gtk_box)
        frame = gtk.Frame(); frame.add(hbox)
        # Add to first column
        vbox_1.pack_start(frame)
        # Store label to allow editing
        self.labels['dtag'] = gtk_value

        # Event Number for Dataset
        gtk_label = gtk.Label('Event')
        gtk_value = gtk.Label('1')
        gtk_box = gtk.EventBox(); gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True); hbox.add(gtk_label); hbox.add(gtk_box)
        frame = gtk.Frame(); frame.add(hbox)
        # Add to first column
        vbox_1.pack_start(frame)
        # Store label to allow editing
        self.labels['e_idx'] = gtk_value

        # Estimated Event Occupancy
        gtk_label = gtk.Label('Est. Occupancy')
        gtk_value = gtk.Label('1.0')
        gtk_box = gtk.EventBox(); gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True); hbox.add(gtk_label); hbox.add(gtk_box)
        frame = gtk.Frame(); frame.add(hbox)
        # Add to first column
        vbox_1.pack_start(frame)
        # Store label to allow editing
        self.labels['e_occ'] = gtk_value

        # Z-Peak for Dataset
        gtk_label = gtk.Label('Z-Peak')
        gtk_value = gtk.Label('0')
        gtk_box = gtk.EventBox(); gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True); hbox.add(gtk_label); hbox.add(gtk_box)
        frame = gtk.Frame(); frame.add(hbox)
        # Add to first column
        vbox_1.pack_start(frame)
        # Store label to allow editing
        self.labels['zpeak'] = gtk_value

        # Z-Peak for Dataset
        gtk_label = gtk.Label('Z-Size')
        gtk_value = gtk.Label('1')
        gtk_box = gtk.EventBox(); gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True); hbox.add(gtk_label); hbox.add(gtk_box)
        frame = gtk.Frame(); frame.add(hbox)
        # Add to first column
        vbox_1.pack_start(frame)
        # Store label to allow editing
        self.labels['csize'] = gtk_value

        return vbox_1

    def _progress_table(self):

        # First Column
        vbox_main = gtk.VBox(spacing=5)
        vbox_main.set_border_width(3)

        # Create title
        title = gtk.Label('Overall Inspection Event/Site Progress:')
        title.set_justify(gtk.JUSTIFY_LEFT)
        frame = gtk.Frame(); frame.add(title)
        # Add to first column
        vbox_main.pack_start(frame)

        # Event Number Name
        gtk_label_1 = gtk.Label('Event')
        gtk_value_1 = gtk.Label(0)
        gtk_label_2 = gtk.Label('of')
        gtk_value_2 = gtk.Label(0)
        # Add values to boxes
        gtk_box_1 = gtk.EventBox(); gtk_box_1.add(gtk_value_1)
        gtk_box_2 = gtk.EventBox(); gtk_box_2.add(gtk_value_2)
        hbox = gtk.HBox(homogeneous=True); hbox.add(gtk_label_1); hbox.add(gtk_box_1); hbox.add(gtk_label_2); hbox.add(gtk_box_2)
        frame = gtk.Frame(); frame.add(hbox)
        # Add to second column
        vbox_main.pack_start(frame)
        # Store label to allow editing
        self.labels['rank_val'] = gtk_value_1
        self.labels['rank_tot'] = gtk_value_2

        # Event Number Name
        gtk_label_1 = gtk.Label('Site')
        gtk_value_1 = gtk.Label(0)
        gtk_label_2 = gtk.Label('of')
        gtk_value_2 = gtk.Label(0)
        # Add values to boxes
        gtk_box_1 = gtk.EventBox(); gtk_box_1.add(gtk_value_1)
        gtk_box_2 = gtk.EventBox(); gtk_box_2.add(gtk_value_2)
        hbox = gtk.HBox(homogeneous=True); hbox.add(gtk_label_1); hbox.add(gtk_box_1); hbox.add(gtk_label_2); hbox.add(gtk_box_2)
        frame = gtk.Frame(); frame.add(hbox)
        # Add to second column
        vbox_main.pack_start(frame)
        # Store label to allow editing
        self.labels['site_val'] = gtk_value_1
        self.labels['site_tot'] = gtk_value_2

        return vbox_main

#############################################################################################
#
# INITIALISE THE MAIN INSPECTOR OBJECT
#
#############################################################################################

work_dir = os.getcwd()
hit_list = os.path.join(work_dir, 'analyses', 'pandda_analyse_events.csv')

inspector = PanddaInspector(csv=hit_list, top_dir=work_dir)
inspector.start_gui()

inspector.refresh_event()
