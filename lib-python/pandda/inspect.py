
import os, sys, glob, copy, time, itertools
import gtk
import pandas

from bamboo.plot import bar
from pandda.constants import PanddaInspectorFilenames, PanddaDatasetFilenames
from pandda.jiffies import pandda_summary
from pandda import inspect_html

def nonmodal_msg(msg):
    """Display an error window - non-modal"""
    d = gtk.MessageDialog(  type    = gtk.MESSAGE_INFO,
                            message_format = msg )
    d.show_all()
    while gtk.events_pending(): gtk.main_iteration(False)
    return d

def modal_msg(msg):
    """Display an error window - model"""
    d = gtk.MessageDialog(  type    = gtk.MESSAGE_INFO,
                            buttons = gtk.BUTTONS_CLOSE,
                            message_format = msg )
    d.run()
    d.destroy()

class PanddaEvent(object):
    def __init__(self, rank, info, top_dir):

        # Key for the pandas table
        self.index = info.name
        # Dataset tag
        self.dtag = info.name[0]
        # Dataset Information
        self.map_resolution  = round(info['analysed_resolution'],2)
        self.map_uncertainty = round(info['map_uncertainty'],2)
        self.rwork_rfree     = (round(info['r_work'],3),round(info['r_free'],3))
        # Event number for the dataset
        self.event_idx = int(info.name[1])
        # Position in the ranked list (1 -> n)
        self.rank = rank
        # Site Number (1 -> m)
        self.site_idx  = int(info['site_idx'])
        # Event Info
        self.est_1_bdc = round(info['1-BDC'], 2)
        # Z statistics
        self.z_peak = round(info['z_peak'], 1)
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
        self.event_map          = os.path.join(self.dataset_dir,    PanddaDatasetFilenames.event_map.format(dtag, self.event_idx, self.est_1_bdc))
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
#            assert os.path.islink(self.fitted_link), 'FILE IS NOT A LINK'
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

    def __init__(self, event_csv, site_csv, top_dir):

        self.log_table  = None
        self.site_table = None

        # List of events from pandda.analyse
        self.site_list = PanddaSiteTracker(parent=self, csv=event_csv, top_dir=top_dir)
        # Handling of coot commands
        self.coot = PanddaMolHandler(parent=self)

        # Working Directory
        self.top_dir = top_dir
        # Output csv from pandda.inspect
        self.output_event_csv = os.path.join(self.top_dir, 'analyses', PanddaInspectorFilenames.event_info)
        self.output_site_csv  = os.path.join(self.top_dir, 'analyses', PanddaInspectorFilenames.site_info)

        # Load previous data or create new table from site list, so we can record inspection data
        self.initialise_output_event_table( in_csv=event_csv, out_csv=self.output_event_csv )
        self.initialise_output_site_table(  in_csv=site_csv,  out_csv=self.output_site_csv  )
        # Save Tables
        self.write_output_csvs()

        # Object to hold the currently loaded event
        self.current_event = None

    def start_gui(self):
        self.gui = PanddaGUI(parent=self)
        self.gui.launch()

    def update_html(self):
        """Update various parts of the class - including output graphs"""

        # Post message in case this takes a long time...
        d = nonmodal_msg('Updating html output...')
        # Plot output graph of site list
        print 'Generating output graphs...'
        plot_vals = self.log_table['z_peak']
        view_vals = self.log_table['Viewed']
        modl_vals = self.log_table['Ligand Placed']
        colr_vals = ['limegreen' if m else 'red' if v else 'blue' for m,v in zip(modl_vals,view_vals)]
        site_idxs = self.log_table['site_idx']
        groups = [list(g[1]) for g in itertools.groupby(range(len(site_idxs)), key=lambda i: site_idxs[i])]
        bar.multiple_bar_plot_over_several_images(
                                f_template  = os.path.join(self.top_dir, 'results_summaries', PanddaInspectorFilenames.inspect_site_graph_mult),
                                plot_vals   = [[plot_vals[i] for i in g] for g in groups],
                                colour_vals = [[colr_vals[i] for i in g] for g in groups]   )
        # Write output html
        print 'Writing output html...'
        inspect_html.write_inspect_html(top_dir=self.top_dir, inspector=self)
        # Destroy update message
        d.destroy()

    def update_gui(self):

        # Update global position
        self.gui.labels['site_val'].set_label(str(self.site_list.site_val))
        self.gui.labels['site_tot'].set_label(str(self.site_list.site_tot))
        self.gui.labels['rank_val'].set_label(str(self.site_list.rank_val))
        self.gui.labels['rank_tot'].set_label(str(self.site_list.rank_tot))
        # Update current event and dataset information
        self.gui.labels['dtag'].set_label(str(self.current_event.dtag))
        self.gui.labels['e_idx'].set_label(str(self.current_event.event_idx))
        self.gui.labels['e_1_bdc'].set_label(str(self.current_event.est_1_bdc))
        self.gui.labels['zpeak'].set_label(str(round(self.current_event.z_peak,3)))
        self.gui.labels['csize'].set_label(str(self.current_event.cluster_size))
        self.gui.labels['map_res'].set_label(str(self.current_event.map_resolution))
        self.gui.labels['map_unc'].set_label(str(self.current_event.map_uncertainty))
        self.gui.labels['rwork_rfree'].set_label('{} / {}'.format(*self.current_event.rwork_rfree))

        # Reset the event comment boxes
        self.gui.objects['event comment text'].set_text(str(self.get_event_log_value('Comment')))

        # Reset the site comment boxes
        self.gui.objects['site name text'].set_text(str(self.site_table.get_value(index=self.current_event.site_idx, col='Name')))
        self.gui.objects['site comment text'].set_text(str(self.site_table.get_value(index=self.current_event.site_idx, col='Comment')))

        # Update the radio buttons - "Interesting"
        if self.get_event_log_value('Interesting') == True:
            self.gui.buttons['tp'].set_active(True)
        else:
            self.gui.buttons['fp'].set_active(True)
        # Update the radio buttons - "Ligand Placed"
        if self.get_event_log_value('Ligand Placed') == True:
            self.gui.buttons['placed'].set_active(True)
        else:
            self.gui.buttons['not placed'].set_active(True)
        # Update the radio buttons - "Ligand Confidence"
        if self.get_event_log_value('Ligand Confidence') == 'High':
            self.gui.buttons['high conf'].set_active(True)
        elif self.get_event_log_value('Ligand Confidence') == 'Medium':
            self.gui.buttons['med conf'].set_active(True)
        else:
            self.gui.buttons['low conf'].set_active(True)

        # Reset the merge button
        self.gui.buttons['merge'].child.set_text("Merge Ligand\nWith Model")

    #-------------------------------------------------------------------------

    def save_current(self):
        self.current_event.write_fitted_model(protein_obj=self.coot.open_mols['p'])
        self.write_output_csvs()

    def reset_current_to_last_model(self):
        close_molecule(self.coot.open_mols['p'])
        if os.path.exists(self.current_event.fitted_link):
            p = read_pdb(self.current_event.fitted_link)
        else:
            p = read_pdb(self.current_event.unfitted_model)
        self.coot.open_mols['p'] = p
        self.write_output_csvs()

    def reset_current_to_orig_model(self):
        close_molecule(self.coot.open_mols['p'])
        p = read_pdb(self.current_event.unfitted_model)
        self.coot.open_mols['p'] = p
        self.write_output_csvs()

    def load_new_event(self, new_event):
        if new_event is None:
            # TODO DO SOMETHING BETTER HERE TODO
            # Instead of loading new event, inform using the same model
            modal_msg(msg='Reloading Same MODEL (Unchanged)')
        else:
            self.current_event = self.coot.load_event(e=new_event)
        self.update_gui()

    def refresh_event(self):
        self.load_new_event(new_event=self.site_list.get_new_current_event())

    def load_next_event(self, skip_unmodelled=False, skip_viewed=False):
        if not self.check_valid_button_click(skip_unmodelled, skip_viewed): return
        new_event = self.site_list.get_next()
        # Loop until we get a modelled one (if requested)
        if (new_event is not None) and skip_unmodelled and (not os.path.exists(new_event.fitted_link)):
            if self.current_event.index == new_event.index:
                modal_msg(msg='No modelled datasets')
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

    def load_dataset(self, dataset_id):
        """Find the next dataset with the given id"""
        new_event = self.site_list.get_next()
        # Check if this is the right dataset
        if (new_event is not None) and (new_event.dtag != dataset_id):
            if self.current_event.index == new_event.index:
                # Check if we've looped around
                modal_msg(msg='No dataset found for this id')
            else:
                # Load the next dataset if this doesn't match
                print 'Event does not match dataset id: {} - {}'.format(dataset_id, new_event.dtag)
                self.load_dataset(dataset_id=dataset_id)
            return
        # Actually load the event
        self.load_new_event(new_event=new_event)

    #-------------------------------------------------------------------------

    def check_valid_button_click(self, skip_unmodelled=None, skip_viewed=None):
        if skip_unmodelled == True:
            pass
        if skip_viewed == True:
            if sum(self.log_table['Viewed']) == len(self.log_table['Viewed']):
                modal_msg(msg='All models have been viewed')
                return False
        return True
    #-------------------------------------------------------------------------

    def initialise_output_event_table(self, in_csv, out_csv):
        """Read in the log table from pandda.analyse and merge with previous results form pandda.inspect"""

        # Read in the log table from pandda_analyse
        self.log_table = pandas.read_csv(in_csv, sep=',', dtype={'dtag':str})
        self.log_table = self.log_table.set_index(['dtag','event_idx'])
        # Create new columns (filled with blanks or defaults)
        self.log_table['Interesting'] = False
        self.log_table['Ligand Placed'] = False
        self.log_table['Ligand Confidence'] = 'Low'
        self.log_table['Comment'] = 'None'
        self.log_table['Viewed'] = False

        if os.path.exists(out_csv):
            print 'Merging with existing pandda_inspect_events.csv...'
            # Output csv already exists from previous run - reload and merge with in_csv
            inspect_prev = pandas.read_csv(out_csv, sep=',', dtype={'dtag':str})
            inspect_prev = inspect_prev.set_index(['dtag','event_idx'])
            # Merge with input table (only on the columns that should be updated)
            self.log_table.update(inspect_prev[['Interesting','Ligand Placed','Ligand Confidence','Comment','Viewed']])

    def initialise_output_site_table(self, in_csv, out_csv):
        """Read in the site definition table from pandda.analyse and merge with previous results form pandda.inspect"""

        # Read in the log table from pandda_analyse
        self.site_table = pandas.read_csv(in_csv, sep=',')
        self.site_table = self.site_table.set_index('site_idx')
        # Create new columns (filled with blanks for the moment)
        self.site_table['Name']     = 'None'
        self.site_table['Comment']  = 'None'

        if os.path.exists(out_csv):
            print 'Merging with existing pandda_inspect_events.csv...'
            # Output csv already exists from previous run - reload and merge with in_csv
            inspect_prev = pandas.read_csv(out_csv, sep=',')
            inspect_prev = inspect_prev.set_index('site_idx')
            # Merge with input table (only on the columns that should be updated)
            self.site_table.update(inspect_prev[['Name','Comment']])

    def get_event_log_value(self, col):
        return self.log_table.get_value(index=self.current_event.index, col=col)

    def set_event_log_value(self, col, value, write=True):
        self.log_table.set_value(index=self.current_event.index, col=col, value=value)
        print '====================>>>'
        print 'EVENT TABLE UPDATED'
        print '====================>>>'
        print 'Current Event:'
        print self.log_table.loc[self.current_event.index]
        print '====================>>>'
        if write: self.write_output_csvs()

    def get_site_log_value(self, col):
        return self.site_table.get_value(index=self.current_event.site_idx, col=col)

    def set_site_log_value(self, col, value, write=True):
        self.site_table.set_value(index=self.current_event.site_idx, col=col, value=value)
        print '====================>>>'
        print 'SITE TABLE UPDATED:'
        print '====================>>>'
        print 'Current Site:'
        print self.site_table.loc[self.current_event.site_idx]
        print '====================>>>'
        if write: self.write_output_csvs()

    def write_output_csvs(self):
        print 'Writing output csv: {}'.format(self.output_event_csv)
        self.log_table.to_csv(self.output_event_csv)
        print 'Writing output csv: {}'.format(self.output_site_csv)
        self.site_table.to_csv(self.output_site_csv)

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

        # Event Map
        o = handle_read_ccp4_map(e.event_map, 0)
        set_last_map_contour_level(2*e.est_1_bdc)
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

            # Set the b-factors of the ligand to 20
            set_b_factor_molecule(l, 20)
            # Set the occupancy of the ligand to 2*(1-bdc)
            all_residue_ids = all_residues(l)
            if all_residue_ids:
                for res_chn, res_num, res_ins in all_residue_ids:
                    set_alt_conf_occ(l, res_chn, res_num, res_ins, [['', 2.0*e.est_1_bdc]])

        return e

class PanddaGUI(object):

    """GUI Class for pandda.inspect"""

    def __init__(self, parent):

        self.parent = parent

        # GUI objects
        self.labels = {}
        self.buttons = {}
        self.objects = {}

    def quit(self):
        self.store()
        self.parent.update_html()
        gtk.main_quit()

    def on_destroy(self,  widget=None, *data):
        self.quit()

    def on_delete(self, widget=None, *data):
        self.quit()

    def launch(self):
        """Launch GUI window"""

        # Create main window
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.connect("delete_event", self.on_delete)
        self.window.connect("destroy_event", self.on_destroy)
        self.window.set_border_width(10)
        self.window.set_default_size(600, 400)
        self.window.set_title("PANDDA inspect")

        # Main VBox object for the window
        main_vbox = gtk.VBox(spacing=5)
        # -----------------------------------------------------
        hbox = gtk.HBox(spacing=5)
        # Create buttones to allow user to quit
        quit_buttons = self._quit_buttons()
        frame = gtk.Frame(); frame.add(quit_buttons)
        hbox.pack_start(frame)
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
        rec_e_buttons = self._record_event_buttons()
        frame = gtk.Frame(); frame.add(rec_e_buttons)
        main_vbox.pack_start(frame)
        # -----------------------------------------------------
        # Create buttones to record meta about the event
        rec_s_buttons = self._record_site_buttons()
        frame = gtk.Frame(); frame.add(rec_s_buttons)
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

        self.buttons['go-to'].connect("clicked", lambda x: [self.store(), self.parent.load_dataset(dataset_id=self.objects['go-to-text'].get_text().strip()), self.objects['go-to-text'].set_text('')])

        # Quit
        self.buttons['quit'].connect("clicked", lambda x: [self.quit()])
        self.buttons['summary'].connect("clicked", lambda x: [self.store(), self.parent.update_html(), os.system('ccp4-python -Qnew -m pandda.jiffies.pandda_summary &')])
#        self.buttons['summary'].connect("clicked", lambda x: [self.store(), self.parent.update_html(), os.system('pandda.show_summary &')])
#        self.buttons['summary'].connect("clicked", lambda x: [self.store(), self.parent.update_html(), pandda_summary.run()])

        self.buttons['updatehtml'].connect("clicked", lambda x: [self.store(), self.parent.update_html()])

        # Structure Buttons
        self.buttons['save'].connect("clicked",  lambda x: self.parent.save_current())
        self.buttons['reload'].connect("clicked", lambda x: self.parent.reset_current_to_last_model())
        self.buttons['reset'].connect("clicked",  lambda x: self.parent.reset_current_to_orig_model())

        # Ligand Buttons
        self.buttons['merge'].connect("clicked", lambda x: [self.buttons['merge'].child.set_text('Already Merged\n(Click to Repeat)'), self.parent.coot.merge_ligand_with_protein()])
        self.buttons['move'].connect("clicked",  lambda x: self.parent.coot.move_ligand_here())

        # Meta Recording buttons
        self.buttons['tp'].connect("clicked",         lambda x: self.parent.set_event_log_value(col='Interesting', value=True))
        self.buttons['fp'].connect("clicked",         lambda x: self.parent.set_event_log_value(col='Interesting', value=False))
        self.buttons['high conf'].connect("clicked",  lambda x: self.parent.set_event_log_value(col='Ligand Confidence', value='High'))
        self.buttons['med conf'].connect("clicked",   lambda x: self.parent.set_event_log_value(col='Ligand Confidence', value='Medium'))
        self.buttons['low conf'].connect("clicked",   lambda x: self.parent.set_event_log_value(col='Ligand Confidence', value='Low'))
        self.buttons['placed'].connect("clicked",     lambda x: self.parent.set_event_log_value(col='Ligand Placed', value=True))
        self.buttons['not placed'].connect("clicked", lambda x: self.parent.set_event_log_value(col='Ligand Placed', value=False))

    def store(self):
        """Record information from the gui to the pandas table in the main object"""
        # Event records
        self.parent.set_event_log_value(col='Comment', value=self.objects['event comment text'].get_text(), write=False)
        self.parent.set_event_log_value(col='Viewed',  value=True, write=False)
        # Site records
        self.parent.set_site_log_value(col='Name',    value=self.objects['site name text'].get_text(),    write=False)
        self.parent.set_site_log_value(col='Comment', value=self.objects['site comment text'].get_text(), write=False)
        # Write csvs only once
        self.parent.write_output_csvs()

    def _navi_buttons_1(self):
        box = gtk.HBox(homogeneous=False, spacing=2)
        box.set_border_width(3)
        # ---
        b = gtk.Button(label="<<< Prev <<<\n(Don't Save Model)")
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        self.buttons['prev'] = b
        box.pack_start(b)
        # ---
        b = gtk.Button(label=">>> Next >>>\n(Don't Save Model)")
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        self.buttons['skip'] = b
        box.pack_start(b)
        # ---
        b = gtk.VSeparator()
        box.pack_start(b, expand=False, padding=5)
        # ---
        b = gtk.Button(label=">>> Next >>>\n(Save Model)")
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        self.buttons['next'] = b
        box.pack_start(b)
        # ---
        return box

    def _navi_buttons_2(self):
        main_box  = gtk.VBox(homogeneous=True, spacing=5)
        # ---
        box1 = gtk.HBox(homogeneous=False, spacing=2)
        box1.set_border_width(3)
        frame = gtk.Frame(); frame.add(box1)
        main_box.pack_start(frame)
        # ---
        box2 = gtk.HBox(homogeneous=True, spacing=2)
        box2.set_border_width(3)
        frame = gtk.Frame(); frame.add(box2)
        main_box.pack_start(frame)
        # ---
        box3 = gtk.HBox(homogeneous=True, spacing=2)
        box3.set_border_width(3)
        frame = gtk.Frame(); frame.add(box3)
        main_box.pack_start(frame)
        # ---
        l = gtk.Label('Go to Dataset:')
        box1.pack_start(l, expand=False, fill=False, padding=5)
        # ---
        e = gtk.Entry(max=200)
        self.objects['go-to-text'] = e
        box1.pack_start(e, expand=True, fill=True, padding=5)
        # ---
        b = gtk.Button(label="Go")
        self.buttons['go-to'] = b
        box1.pack_start(b, expand=False, fill=False, padding=5)
        # ---
        b = gtk.Button(label="<<< Go to Prev Site <<<")
        self.buttons['prev-site'] = b
        box2.add(b)
        # ---
        b = gtk.Button(label=">>> Go to Next Site >>>")
        self.buttons['next-site'] = b
        box2.add(b)
        # ---
        b = gtk.Button(label=">>> Go to Next Unviewed >>>")
        self.buttons['next-unviewed'] = b
        box3.add(b)
        # ---
        b = gtk.Button(label=">>> Go to Next Modelled >>>")
        self.buttons['next-modelled'] = b
        box3.add(b)
        # ---
        return main_box

    def _ligand_buttons(self):
        box1 = gtk.VBox(homogeneous=True, spacing=2)
        box1.set_border_width(3)
        # ---
        b = gtk.Button(label="Merge Ligand\nWith Model")
        b.child.set_line_wrap(True)
        b.child.props.width_chars = 15
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        self.buttons['merge'] = b
        box1.add(b)
        # ---
        b = gtk.Button(label="Move New Ligand Here")
        b.child.set_line_wrap(True)
        b.child.props.width_chars = 10
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        self.buttons['move'] = b
        box1.add(b)
        # ---
        box2 = gtk.VBox(homogeneous=True, spacing=2)
        box2.set_border_width(3)
        # ---
        b = gtk.Button(label="Save Model")
        b.child.set_line_wrap(True)
        b.child.props.width_chars = 10
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        self.buttons['save'] = b
        box2.add(b)
        # ---
        b = gtk.Button(label="Reload Last Saved Model")
        b.child.set_line_wrap(True)
        b.child.props.width_chars = 15
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        self.buttons['reload'] = b
        box2.add(b)
        # ---
        b = gtk.Button(label="Reset to Unfitted Model")
        b.child.set_line_wrap(True)
        b.child.props.width_chars = 15
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        self.buttons['reset'] = b
        box2.add(b)
        # ---
        hbox_main = gtk.HBox(spacing=5)
        hbox_main.pack_start(box1)
        hbox_main.pack_start(gtk.VSeparator(), expand=False, padding=5)
        hbox_main.pack_start(box2)

        return hbox_main

    def _record_event_buttons(self):
        # ---------------------------------------------
        hbox_1 = gtk.HBox(homogeneous=False, spacing=5)
        # ---
        hbox_1.pack_start(gtk.HBox(), expand=False, fill=False, padding=10)
        # ---
        l = gtk.Label('Event Comment:')
        hbox_1.pack_start(l, expand=False, fill=False, padding=5)
        # ---
        e = gtk.Entry(max=200)
        self.objects['event comment text'] = e
        hbox_1.pack_start(e, expand=True, fill=True, padding=5)
        # ---
        hbox_1.pack_start(gtk.HBox(), expand=False, fill=False, padding=10)
        # ---------------------------------------------
        hbox_2 = gtk.HBox(homogeneous=True, spacing=5)
        hbox_2.set_border_width(3)
        vbox_1_1 = gtk.VBox(homogeneous=True, spacing=2)
        vbox_1_2 = gtk.VBox(homogeneous=True, spacing=2)
        vbox_1_3 = gtk.VBox(homogeneous=True, spacing=2)
        hbox_2.add(vbox_1_1); hbox_2.add(vbox_1_2); hbox_2.add(vbox_1_3)
        # ---
        b = gtk.RadioButton(label="Mark Event as Interesting")
        self.buttons['tp'] = b
        vbox_1_1.add(b)
        # ---
        b = gtk.RadioButton(label="Mark Event as Not Interesting", group=b)
        self.buttons['fp'] = b
        vbox_1_1.add(b)
        # ---
        b = gtk.RadioButton(label="Ligand Placed")
        self.buttons['placed'] = b
        vbox_1_2.add(b)
        # ---
        b = gtk.RadioButton(label="No Ligand Placed", group=b)
        self.buttons['not placed'] = b
        vbox_1_2.add(b)
        # ---
        b = gtk.RadioButton(label="Model: High Confidence")
        self.buttons['high conf'] = b
        vbox_1_3.add(b)
        # ---
        b = gtk.RadioButton(label="Model: Medium Confidence", group=b)
        self.buttons['med conf'] = b
        vbox_1_3.add(b)
        # ---
        b = gtk.RadioButton(label="Model: Low Confidence", group=b)
        self.buttons['low conf'] = b
        vbox_1_3.add(b)
        # ---------------------------------------------
        vbox_main = gtk.VBox(spacing=0)
        vbox_main.pack_start(gtk.Label('Record Event Information (this event only)'), expand=False, fill=False, padding=5)
        vbox_main.pack_start(hbox_1, padding=0)
        vbox_main.pack_start(hbox_2, padding=5)

        return vbox_main

    def _record_site_buttons(self):
        # ---------------------------------------------
        hbox_1 = gtk.HBox(homogeneous=False, spacing=5)
        # ---
        hbox_1.pack_start(gtk.HBox(), expand=False, fill=False, padding=10)
        # ---
        l = gtk.Label('Name:')
        l.set_width_chars(10)
        hbox_1.pack_start(l, expand=False, fill=False, padding=5)
        # ---
        e = gtk.Entry(max=200)
        self.objects['site name text'] = e
        hbox_1.pack_start(e, expand=True, fill=True, padding=5)
        # ---
        hbox_1.pack_start(gtk.HBox(), expand=False, fill=False, padding=10)
        # ---------------------------------------------
        hbox_2 = gtk.HBox(homogeneous=False, spacing=5)
        # ---
        hbox_2.pack_start(gtk.HBox(), expand=False, fill=False, padding=10)
        # ---
        l = gtk.Label('Comment:')
        l.set_width_chars(10)
        hbox_2.pack_start(l, expand=False, fill=False, padding=5)
        # ---
        e = gtk.Entry(max=200)
        self.objects['site comment text'] = e
        hbox_2.pack_start(e, expand=True, fill=True, padding=5)
        # ---
        hbox_2.pack_start(gtk.HBox(), expand=False, fill=False, padding=10)
        # ---------------------------------------------
        vbox_main = gtk.VBox(spacing=0)
        vbox_main.pack_start(gtk.Label('Record Site Information (for all events with this site)'), expand=False, fill=False, padding=5)
        vbox_main.pack_start(hbox_1, padding=5)
        vbox_main.pack_start(hbox_2, padding=5)

        return vbox_main

    def _quit_buttons(self):
        # ---------------------------------------------
        vbox_1 = gtk.VBox(spacing=5)
        vbox_1.set_border_width(3)
        # ---
        b = gtk.Button(label="Quit")
        self.buttons['quit'] = b
        vbox_1.pack_start(b)
        # ---
        vbox_1.pack_start(gtk.HSeparator(), expand=False, padding=2)
        # ---
        b = gtk.Button(label="Summary")
        self.buttons['summary'] = b
        vbox_1.pack_start(b)
        # ---
        vbox_1.pack_start(gtk.HSeparator(), expand=False, padding=2)
        # ---
        b = gtk.Button(label="Update HTML")
        self.buttons['updatehtml'] = b
        vbox_1.pack_start(b)
        # ---
        hbox_main = gtk.HBox(spacing=5)
        hbox_main.pack_start(vbox_1)

        return hbox_main

    def _event_info_table(self):

        # Main box
        vbox_main = gtk.VBox()
        vbox_main.set_border_width(3)
        #  Pack sub-boxes
        hbox_sub_1 = gtk.HBox()
        hbox_sub_2 = gtk.HBox()
        vbox_main.pack_start(hbox_sub_1)
        vbox_main.pack_start(hbox_sub_2)

        ##############
        # HBOX SUB 1 #
        ##############

        # Dataset Name
        gtk_label = gtk.Label('Dataset ID')
        gtk_value = gtk.Label('None')
        gtk_box = gtk.EventBox(); gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True); hbox.add(gtk_label); hbox.add(gtk_box)
        frame = gtk.Frame(); frame.add(hbox)
        # Add to first box
        hbox_sub_1.pack_start(frame)
        # Store label to allow editing
        self.labels['dtag'] = gtk_value

        ##############
        # HBOX SUB 2 #
        ##############

        vbox_1 = gtk.VBox()
        vbox_2 = gtk.VBox()
        hbox_sub_2.pack_start(vbox_1)
        hbox_sub_2.pack_start(vbox_2)

        ##########
        # VBOX 1 #
        ##########

        # Create title
        title = gtk.Label('Event Information:')
        title.set_justify(gtk.JUSTIFY_LEFT)
        frame = gtk.Frame(); frame.add(title)
        # Add to first column
        vbox_1.pack_start(frame)

        # Event Number for Dataset
        gtk_label = gtk.Label('Event #')
        gtk_value = gtk.Label('-')
        gtk_box = gtk.EventBox(); gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True); hbox.add(gtk_label); hbox.add(gtk_box); hbox.set_border_width(3)
        frame = gtk.Frame(); frame.add(hbox)
        # Add to first column
        vbox_1.pack_start(frame)
        # Store label to allow editing
        self.labels['e_idx'] = gtk_value

        # Estimated Event Background Correction
        gtk_label = gtk.Label('1 - BDC')
        gtk_value = gtk.Label('-')
        gtk_box = gtk.EventBox(); gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True); hbox.add(gtk_label); hbox.add(gtk_box); hbox.set_border_width(3)
        frame = gtk.Frame(); frame.add(hbox)
        # Add to first column
        vbox_1.pack_start(frame)
        # Store label to allow editing
        self.labels['e_1_bdc'] = gtk_value

        # Z-Peak for Dataset
        gtk_label = gtk.Label('Z-blob Peak')
        gtk_value = gtk.Label('-')
        gtk_box = gtk.EventBox(); gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True); hbox.add(gtk_label); hbox.add(gtk_box); hbox.set_border_width(3)
        frame = gtk.Frame(); frame.add(hbox)
        # Add to first column
        vbox_1.pack_start(frame)
        # Store label to allow editing
        self.labels['zpeak'] = gtk_value

        # Z-Peak for Dataset
        gtk_label = gtk.Label('Z-blob Size')
        gtk_value = gtk.Label('-')
        gtk_box = gtk.EventBox(); gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True); hbox.add(gtk_label); hbox.add(gtk_box); hbox.set_border_width(3)
        frame = gtk.Frame(); frame.add(hbox)
        # Add to first column
        vbox_1.pack_start(frame)
        # Store label to allow editing
        self.labels['csize'] = gtk_value

        ##########
        # VBOX 2 #
        ##########

        # Create title
        title = gtk.Label('Dataset Information:')
        title.set_justify(gtk.JUSTIFY_LEFT)
        frame = gtk.Frame(); frame.add(title)
        # Add to second column
        vbox_2.pack_start(frame)

        # Resolution for Dataset
        gtk_label = gtk.Label('Resolution')
        gtk_value = gtk.Label('-')
        gtk_box = gtk.EventBox(); gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True); hbox.add(gtk_label); hbox.add(gtk_box); hbox.set_border_width(3)
        frame = gtk.Frame(); frame.add(hbox)
        # Add to second column
        vbox_2.pack_start(frame)
        # Store label to allow editing
        self.labels['map_res'] = gtk_value

        # Map Uncertainty for Dataset
        gtk_label = gtk.Label('Map Uncertainty')
        gtk_value = gtk.Label('-')
        gtk_box = gtk.EventBox(); gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True); hbox.add(gtk_label); hbox.add(gtk_box); hbox.set_border_width(3)
        frame = gtk.Frame(); frame.add(hbox)
        # Add to second column
        vbox_2.pack_start(frame)
        # Store label to allow editing
        self.labels['map_unc'] = gtk_value

        # R-Free/R-Work for Dataset
        gtk_label = gtk.Label('R-Free / R-Work')
        gtk_value = gtk.Label('-')
        gtk_box = gtk.EventBox(); gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True); hbox.add(gtk_label); hbox.add(gtk_box); hbox.set_border_width(3)
        frame = gtk.Frame(); frame.add(hbox)
        # Add to second column
        vbox_2.pack_start(frame)
        # Store label to allow editing
        self.labels['rwork_rfree'] = gtk_value

        # Currently Blank
        gtk_label = gtk.Label('-')
        gtk_value = gtk.Label('-')
        gtk_box = gtk.EventBox(); gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True); hbox.add(gtk_label); hbox.add(gtk_box); hbox.set_border_width(3)
        frame = gtk.Frame(); frame.add(hbox)
        # Add to second column
        vbox_2.pack_start(frame)
        # Store label to allow editing
        self.labels['blank'] = gtk_value

        return vbox_main

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

if __name__=='__main__':

#    import os
#    from pandda.inspect_main import PanddaInspector
    from pandda.constants import PanddaAnalyserFilenames

    #############################################################################################
    #
    # CHANGE COOT SETTINGS
    #
    #############################################################################################

    try:
        set_nomenclature_errors_on_read("ignore")
        set_colour_map_rotation_for_map(0)
        set_colour_map_rotation_on_read_pdb(0)
        set_recentre_on_read_pdb(0)
    except:
        pass

    #############################################################################################
    #
    # FIND THE INPUT CSV FILES
    #
    #############################################################################################

    work_dir = os.getcwd()
    hit_list = os.path.join(work_dir, 'analyses', PanddaAnalyserFilenames.event_info)
    site_csv = os.path.join(work_dir, 'analyses', PanddaAnalyserFilenames.site_info)

    #############################################################################################
    #
    # INITIALISE THE MAIN INSPECTOR OBJECT
    #
    #############################################################################################

    inspector = PanddaInspector(event_csv=hit_list, site_csv=site_csv, top_dir=work_dir)
    inspector.start_gui()
    inspector.refresh_event()
