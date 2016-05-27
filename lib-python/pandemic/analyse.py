from __future__ import print_function

import os, sys, glob, time, gc

from scipy import * # So that it doesn't crash...

import numpy

from scitbx.array_family import flex

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

    # Create blank tables and masks
    pandemic.initialise_dataset_masks_and_tables(pandda=pandda)
    pandemic.initialise_residue_masks_and_tables(pandda=pandda)
    # Populate the basic fields
#    pandemic.populate_residue_masks_and_tables(pandda=pandda)

# ============================================================================>
#
###                 PanDEMIC Processing Functions
#
# ============================================================================>

def pandemic_main_loop(pandemic, pandda):

    # Mask each residue to find the block to analyse
    residue_grid_hash = pandemic.find_residue_locales(pandda=pandda)

    pandda.settings.cpus=7
    pandemic.params.resolution = 1.65

    # Load and analyse the variation in the maps
    map_holder_list = pandemic.select_datasets_and_load_maps( pandda = pandda,
                                                              high_res_large_cutoff = pandemic.params.resolution,
                                                              high_res_small_cutoff = 0 )

    all_pcas = {}

    for i_res, res in enumerate(pandda.reference_dataset().calpha_labels()):

        print('\rExtracting map information for residue {}'.format(''.join(res)), end=''); sys.stdout.flush()

        # Get the indices for this residue
        res_idxs = flex.size_t(residue_grid_hash[res])

        # Output array of all of the map values
        map_value_array = numpy.empty((len(map_holder_list.all()), len(res_idxs)))

        # Iterate through the maps
        for i_mh, mh in enumerate(map_holder_list.all()):

            # Get the dataset handler
            dh = mh.parent

            print('\rExtracting map information for residue {} from dataset {}          '.format(res, dh.tag), end=''); sys.stdout.flush()

            # Select the map values for this residue block
            map_values = mh.map.select(res_idxs)

            # Put the values in the complete array
            map_value_array[i_mh,:] = map_values.as_numpy_array()

#            print('')
#            print(i_mh, dh.tag, map_values.as_numpy_array())
#            print('')

            ########################################################
            # Write out maps for each dataset?
            ########################################################
            if i_res==0:
                from pandda.misc import write_array_to_map
                write_array_to_map( output_file=pandemic.output_handler.get_file(file_tag='dataset_map').format(dh.tag),
                                    map_data=mh.map,
                                    grid=pandda.reference_grid() )
            ########################################################

            ########################################################
            # Write out maps for each residue?
            ########################################################
            if i_mh==0:
                b = flex.bool(mh.map.size()).set_selected(res_idxs, True)
                res_map = flex.double(b.accessor()).set_selected(b, map_values)
                res_map.reshape(mh.map.accessor())
                from pandda.misc import write_array_to_map
                write_array_to_map( output_file=pandemic.output_handler.get_file(file_tag='residue_map').format(dh.tag, ''.join(res)),
                                    map_data=res_map,
                                    grid=pandda.reference_grid() )
            ########################################################

        def perform_pandemic_pca(array):
            """Perform PCA on a map value array to determine the heterogeneity. Each row is a different observation"""

            from sklearn.decomposition import PCA

            pca = PCA()
            pca.fit(array)

            return pca

        # Perform PCA to calculate the number of conformers
        results = perform_pandemic_pca(map_value_array)

        all_pcas[res] = results

        pandemic.log('\rExtracting map information for residue {}: Done                              '.format(''.join(res)))

    pandemic.log('Residue-by-residue electron densities across the datasets have been extracted')

    confs = {}

    for rg in pandda.reference_dataset().hierarchy().residue_groups():
        print(rg)
        confs[(rg.parent().id, rg.resid())] = len(rg.conformers())


    for res in pandda.reference_dataset().calpha_labels():
        pca = all_pcas[res]

        print('====================================>>>')
        print(res)
        print('====================================>>>')
        print('Number of Conformers: {}'.format(confs[res]))
        print('====================================>>>')
        print((pca.explained_variance_ratio_/pca.explained_variance_ratio_[0]).round(3))
        print('====================================>>>')
        from ascii_graph import Pyasciigraph
        g=Pyasciigraph()
        graph_data = [(i+1, val) for i,val in enumerate((pca.explained_variance_ratio_/pca.explained_variance_ratio_[0]).round(3))][0:10]
        for l in g.graph(label='Sorted PCA components (Ascending Order)', data=graph_data, sort=0):
            pandemic.log(l.replace(u"\u2588", '=').replace('= ','> ') True)

    from IPython import embed; embed(); raise SystemExit()



    # Variation of the centre of mass of the protein -- measure of non-isomorphism


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
            print('<<< Pandemic exited before being initialised >>>')
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

