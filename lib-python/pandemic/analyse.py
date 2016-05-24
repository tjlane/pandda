import os, sys, glob, time, gc

from giant.jiffies import extract_params_default

from pandda.analyse_main import PanddaMultiDatasetAnalyser

from pandemic import welcome
from pandemic.phil import pandemic_phil
from pandemic.analyse_main import PandemicMultiDatasetAnalyser

# ============================================================================>
#
###                 PanDDA Initialisation Functions
#
# ============================================================================>

def sanitise_params(params):
    """Ensure continuity in the params as some are redundant"""

    # Change the pandda parameters to match the pandemic parameters
    params.pandda.output.out_dir = params.pandemic.input.pandda_dir

    # This must be set to True so that "old" datasets are analysed
    params.pandda.method.reprocess_existing_datasets = True

def load_pandda_for_pandemic(params, pandemic):
    """Load a pre-calculated pandda using the parameters for pandemic"""


    pandemic.log('===================================>>>', True)
    pandemic.log('==>>>', True)
    pandemic.log('==>>>   Loading the existing PanDDA', True)
    pandemic.log('==>>>', True)
    pandemic.log('===================================>>>', True)

    # Initialise the pandda object
    pandda = PanddaMultiDatasetAnalyser(params)
    # Transfer custom objects from pandemic > pandda
    pandda.log = pandemic.log
    # Load the elements that we need
    pandda.load_pickled_objects()
    # ===============================================================================>
    # CHANGE INTO OUTPUT DIRECTORY (NEEDED FOR READING PICKLES)
    # ===============================================================================>
    os.chdir(pandda.out_dir)
    # Rebuild the pandda masks (quick)
    pandda.initialise_dataset_masks_and_tables()
    pandda.check_loaded_datasets(datasets=pandda.datasets.all())
    pandda.filter_datasets_1()
    pandda.collate_dataset_variables()
    pandda.filter_datasets_2()
    pandda.calculate_dataset_rmsds_to_reference()

    return pandda

# ============================================================================>
#
###                 PanDEMIC Initialisation Functions
#
# ============================================================================>

def pandemic_setup(pandemic, pandda):
    """Prepare for analysis"""

    pandemic.initialise_dataset_masks_and_tables(pandda)
    pandemic.initialise_residue_masks_and_tables(pandda)

# ============================================================================>
#
###                 PanDEMIC Processing Functions
#
# ============================================================================>

def pandemic_main_loop(pandemic, pandda):

    # Variation of the centre of mass of the protein -- measure of non-isomorphism

    # Load and analyse the variation in the maps
    map_holder_list = pandemic.select_datasets_and_load_maps( pandda = pandda,
                                                              high_res_large_cutoff = pandemic.params.resolution,
                                                              high_res_small_cutoff = 0 )

# ============================================================================>
#
###                 PanDEMIC Output + Wrap-Up Functions
#
# ============================================================================>

def pandemic_end(pandemic, pandda):
    pass

# ============================================================================>
#
###                 PanDEMIC Command-Line Function
#
# ============================================================================>

def pandemic_analyse_main(args):
    """Run the PanDEMIC algorithm on the output of a previous pandda analysis"""

    working_phil = extract_params_default(master_phil=pandemic_phil, args=args)
    working_params = working_phil.extract()

    sanitise_params(params=working_params)

    try:
        # ============================================================================>
        #####
        # Initialise a new pandemic instance
        #####
        # ============================================================================>
        pandemic = PandemicMultiDatasetAnalyser(params=working_params)
        pandemic.run_analysis_init()
        # ============================================================================>
        #####
        # Initialise a pre-existing pandda analysis
        #####
        # ============================================================================>
        # Initialise the main object & load the pickled datasets
        pandda = load_pandda_for_pandemic(params=working_params, pandemic=pandemic)
        # ============================================================================>
        #####
        # Use the PanDDA object to populate a new PanDEMIC object
        #####
        # ============================================================================>
        pandemic_setup(pandemic=pandemic, pandda=pandda)
        # ============================================================================>
        #####
        # Run the main analysis loop - now of PanDEMIC
        #####
        # ============================================================================>
        pandemic_main_loop(pandemic=pandemic, pandda=pandda)
        # ============================================================================>
        #####
        # Write summaries and post-process
        #####
        # ============================================================================>
        pandemic_end(pandemic=pandemic, pandda=pandda)
        # ============================================================================>
        #####
        # End
        #####
        # ============================================================================>
    except KeyboardInterrupt:
        raise
    except AssertionError: # TODO REMOVE THIS BEFORE DSITRIBUTION
        raise
    except SystemExit:
        try:
            pandemic.exit(error=False)
        except:
            print '<<< Pandemic exited before being initialised >>>'
    except:
        pandemic.exit(error=True)
    else:
        pandemic.exit(error=False)

    return pandemic

# ============================================================================>
#
#   COMMAND LINE RUN
#
# ============================================================================>

if __name__ == '__main__':

    welcome()
    pandemic = pandemic_analyse_main(args=sys.argv[1:])

