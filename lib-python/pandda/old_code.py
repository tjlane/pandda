
def run_custom_analyses(pandda):
    """If a list of hits is available, test to see whether the pandda identified them"""

    try:

        # ============================================================================>
        # ======================================>
        # Manual Analyses - These only need processing once
        # ======================================>
        # ============================================================================>
        analyses_dir = pandda.output_handler.get_dir('analyses')
        # ============================================================================>

        if 0:
            pandda.log('===================================>>>')
            pandda.log('Calculating Deviations of C-alphas between structures')

            rms = lambda vals: numpy.sqrt(numpy.mean(numpy.abs(vals)**2))
            norm = lambda vals: numpy.sqrt(numpy.sum(numpy.abs(vals)**2))

            # Pull all c-alpha sites for each structure
            all_sites = numpy.array([d.transform_points_to_reference(d.get_calpha_sites()) for d in pandda.datasets.all()])
            # Calculate the mean x,y,z for each c-alpha
            mean_sites = numpy.mean(all_sites, axis=0)
            # Differences from the mean
            diff_sites = all_sites - mean_sites
            # Euclidean norms of the distances moved
            diff_norms = numpy.apply_along_axis(norm, axis=2, arr=diff_sites)

            with open(os.path.join(analyses_dir,'calpha_variation.csv'), 'w') as fh:
                for row in diff_norms:
                    out_list = row.round(3).tolist()
                    out_line = ', '.join(map(str,out_list)) + '\n'
                    fh.write(out_line)

            pandda.log('Largest deviation from the mean site: {!s}'.format(diff_norms.max()))
            pandda.log('Average deviation from the mean site: {!s}'.format(diff_norms.mean()))

        # ============================================================================>

        if 0:
            pandda.log('===================================>>>')
            pandda.log('Clustering the Refined Structures')

            distance_matrix = []
            for d1 in pandda.datasets.all():
               distance_matrix.append([d1.transform_points_to_reference(d1.get_calpha_sites()).rms_difference(d2.transform_points_to_reference(d2.get_calpha_sites())) for d2 in pandda.datasets.all()])

            distance_matrix = numpy.array(distance_matrix)

            with open(os.path.join(analyses_dir,'calpha_distance_matrix.csv'), 'w') as fh:
                for row in distance_matrix:
                    out_list = row.round(3).tolist()
                    out_line = ', '.join(map(str,out_list)) + '\n'
                    fh.write(out_line)

        # ============================================================================>

    except:
        pandda.log('FAILURE DURING CUSTOM ANALYSES', True)
        raise

    return 0

# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# TODO REIMPLEMENT THIS IF IT CAN BE DONE SO THAT IT LOOKS NICE
#
#                # Image Blobs in the datasets using ccp4mg
#                if d_handler.hit_clusters:
#
#                    pandda.log('===================================>>>', True)
#                    pandda.log('Imaging blobs in Dataset {!s}'.format(d_handler.tag), True)
#
#                    sorted_blob_indices = d_handler.hit_clusters.sort(sorting_function=max)
#                    blob_num = len(sorted_blob_indices)
#
#                    for blob_rank, blob_grid_peak in enumerate(d_handler.hit_clusters.get_centroids(indices=sorted_blob_indices)):
#
#                        # Only produce a certain number of images
#                        if blob_rank == pandda.params.blob_search.blobs_to_image:
#                            break
#
#                        status_bar(n=blob_rank, n_max=blob_num)
#
#                        blob_cart_peak = [g*pandda.reference_grid().grid_spacing() for g in blob_grid_peak]
#
#                        # Make images of the blob
#                        pandda.image_blob(  script    = d_handler.output_handler.get_file('ccp4mg_script'),
#                                            image     = d_handler.output_handler.get_file('ccp4mg_png'),
#                                            d_handler = d_handler,
#                                            point_no  = blob_rank+1,
#                                            point     = blob_cart_peak,
#                                            towards   = centre_of_mass
#                                        )
#
#                    status_bar(n=blob_num, n_max=blob_num)
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

