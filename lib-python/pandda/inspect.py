import os

from pandda.inspect_main import PanddaInspector
from pandda.constants import PanddaAnalyserFilenames

if __name__=='__main__':

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
