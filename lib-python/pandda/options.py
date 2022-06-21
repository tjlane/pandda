import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys


class Options(object):

    def __init__(self, config):

        self.multiprocessing_objects(config=config)
        self.plotting_objects(config=config)
        self.preprocess_objects(config=config)
        self.process_objects(config=config)
        self.postprocess_objects(config=config)
        self.output_objects(config=config)

    def multiprocessing_objects(self, config):

        from giant import processors

        # Parameters for dispatched jobs
        self.remote_nodes = (
            config.processing.remote_nodes
            )
        self.cpus_per_remote_node = (
            config.processing.cpus_per_remote_node
            )
        self.h_vmem = (
            config.processing.h_vmem
            )
        self.m_mem_free = (
            config.processing.m_mem_free
            )

        # Parameters for local jobs
        self.local_cpus = (
            config.processing.cpus
            )

        # Numebr of cpus to use within shell processing
        process_shells_n_cpus = None

        #####
        # How to process the shells (outer loop)

        if config.processing.process_shells == "serial":

            self.process_shells = processors.Processor()
            process_shells_n_cpus = self.local_cpus

        elif config.processing.process_shells == "luigi":

            self.process_shells = processors.ProcessorLuigi(
                jobs = self.remote_workers,
                parallel_env = "smp",
                run_locally = False,
                n_cpu = self.remote_cpus_per_worker,
                h_vmem = self.h_vmem,
                m_mem_free = self.m_mem_free,
            )
            process_shells_n_cpus = self.remote_cpus_per_worker

        else:

            raise NotImplemented()

        #####

        assert (self.process_shells is not None)

        assert (process_shells_n_cpus is not None)

        #####
        # Processor for within each shell (inner loop)

        if config.processing.backend == "serial":

            self.process_in_shell = processors.Processor()

        elif config.processing.backend == "parallel_joblib":

            self.process_in_shell = processors.ProcessorJoblib(
                n_cpus = process_shells_n_cpus,
                )

        else:

            raise NotImplemented()

        #####

        assert (self.process_in_shell is not None)

        #####
        # Local jobs (e.g. alignment)

        if config.processing.backend == "serial":

            self.process_general = processors.Processor()

        elif config.processing.backend == "parallel_joblib":

            self.process_general = processors.ProcessorJoblib(
                n_cpus = self.local_cpus,
                backend = 'multiprocessing', # loky "broken" on mac
                )

        else:

            raise NotImplemented()

        #
        #####

        assert (self.process_general is not None)

        #####

    def plotting_objects(self, config):

        import pandda.graphs
        from pandda.graphs import configure
        plot_config = configure()
        assert pandda.graphs.config is plot_config

    def preprocess_objects(self, config):

        #####
        # Dataloader
        from giant.mulch.dataloader import MultiDatasetDataloader
        dataloader = MultiDatasetDataloader(
            data_dirs = config.input.data_dirs,
            pdb_style = config.input.pdb_style,
            mtz_style = config.input.mtz_style,
            pdb_regex = config.input.regex.pdb_regex,
            mtz_regex = config.input.regex.mtz_regex,
            dir_regex = config.input.regex.dir_regex,
            only_datasets = config.input.flags.only_datasets, # !!!
            ignore_datasets = config.input.flags.ignore_datasets, # !!!
            dataset_prefix = config.output.dataset_prefix,
            # lig_style = config.input.lig_style, # !!!
        )

        #####
        # Dataset partitioner
        from giant.mulch.partitioners import TestTrainPartitioner
        from giant.mulch.selectors import (
            SortedDatasetSelector,
            #RandomDatasetSelector,
            SoulmateSelector,
            )
        from giant.mulch.sorters import (
            HighResolutionSorter,
            )

        self.partition_test_train = TestTrainPartitioner(
            test = config.input.flags.test,
            train = config.input.flags.train,
            not_test = config.input.flags.not_test,
            not_train = config.input.flags.not_train,
            train_selector = SortedDatasetSelector(
                max_datasets = config.params.statistical_maps.max_build_datasets,
                sort_datasets = HighResolutionSorter(),
                ),
            # train_selector = RandomDatasetSelector(
            #     max_datasets = config.params.statistical_maps.max_build_datasets,
            #     ),
            test_selector = (
                SoulmateSelector( # Finds "the one" and won't select anything else
                    sort_datasets = HighResolutionSorter(),
                    )
                if (config.output.output_maps_for == 'first_dataset_only')
                else None
                ),
            )

        #####
        # Initial dataset filter
        from giant.mulch.filters import (
            DatasetFilterGroup, 
            RValueFilter,
            HighResolutionFilter,
            )
        initial_filter = DatasetFilterGroup(
            filters = [
                # Remove terrible datasets
                RValueFilter(
                    max_rfree = config.params.filtering.max_rfree,
                    # max_rwork = !!!,
                    ),
                # This filter saves loading datasets that would not be analysed anyway!
                HighResolutionFilter(
                    high_resolution_cutoff = config.params.analysis.high_res_lower_limit,
                    ),
                ],
            )

        # Filter by input structure
        if (config.input.filter.pdb is not None):
            from giant.dataset import ModelAndData
            filter_dataset = ModelAndData.from_file(
                structure_filename = config.input.filter.pdb,
                )
            from giant.mulch.filters import DefaultDatasetFilter
            initial_filter.filters.append(
                DefaultDatasetFilter(
                    same_space_group_only = True,
                    similar_models_only = config.params.filtering.similar_models_only,
                    max_rfree = None,
                    max_rwork = None,
                    reference_dataset = filter_dataset,
                    )
                )

        #####
        # Select reference
        from giant.mulch.reference import DefaultReferenceSelector
        select_reference = DefaultReferenceSelector(
            pdb_path = config.input.reference.pdb,
            mtz_path = config.input.reference.mtz,
        )

        #####
        # Filter with reference
        from giant.mulch.filters import DefaultDatasetFilter
        filter_with_reference = DefaultDatasetFilter(
            same_space_group_only = config.params.filtering.same_space_group_only,
            similar_models_only = config.params.filtering.similar_models_only,
            # max_rfree = config.params.filtering.max_rfree,
            # max_rwork = !!!,
        )

        #####
        # Input data checks
        from giant.mulch.checks import MultiDatasetChecker, DiffractionArrayChecker
        check_input_data = MultiDatasetChecker(
            dataset_checker = DiffractionArrayChecker(
                possible_structure_factor_pairs = config.params.diffraction_data.structure_factors,
                required_structure_factor_pairs = None,
                check_completeness_until = config.params.diffraction_data.checks.low_resolution_completeness,
                check_for_invalid_values = config.params.diffraction_data.checks.all_data_are_valid_values,
            )
        )

        #####
        # Get input miller array
        from giant.mulch.xray.get_data import LoadFirstValidMillerArrayFromListOfColumnOptions
        get_input_miller_array = LoadFirstValidMillerArrayFromListOfColumnOptions(
            structure_factor_pairs = config.params.diffraction_data.structure_factors,
            )

        #####
        # Scale miller arrays
        if config.params.diffraction_data.scaling.apply_b_factor_scaling:
            from giant.mulch.xray.scale_data import (
                GetMultiDatasetMillerArrayScalers,
                GetIsotropicMillerArrayScaler,
                )
            get_dataset_scaler = GetMultiDatasetMillerArrayScalers(
                get_miller_array = get_input_miller_array,
                get_scaler = GetIsotropicMillerArrayScaler(),
                )
        else:
            get_dataset_scaler = None

        #####
        # Align datasets
        from giant.mulch.transform.align import AlignDatasets
        align_datasets = AlignDatasets(
            method = config.params.alignment.method,
            processor = self.process_general,
        )

        #####
        # Extract dataset statistics
        from pandda.preprocess.statistics import PanddaExtractDatasetStatistics
        extract_dataset_statistics = PanddaExtractDatasetStatistics(
            max_scaling_z_score = config.params.excluding.max_wilson_plot_z_score,
            max_wilson_z_score = config.params.excluding.max_wilson_plot_z_score,
            )

        #####
        #
        # Output Objects
        #
        from pandda.preprocess.copy_input_files import CopyPanddaInputFiles
        copy_input_files = CopyPanddaInputFiles(
            output_dir = (
                config.output.out_dir / "processed_datasets"
                ),
            )
        from pandda.preprocess.output import WritePanddaDatasetSummary
        write_initial_dataset_summary = WritePanddaDatasetSummary(
            output_dir = (
                config.output.out_dir
                ),
            dataset_dir = (
                config.output.out_dir / "processed_datasets"
                ),
            )

        #####
        #
        # Combine into single task
        #
        from pandda.preprocess import DatasetInitialiser
        self.dataset_initialiser = DatasetInitialiser(
            dataloader = dataloader,
            partitioner = self.partition_test_train,
            initial_filter = initial_filter,
            select_reference = select_reference,
            filter_with_reference = filter_with_reference,
            check_input_data = check_input_data,
            calculate_scalings = get_dataset_scaler,
            get_input_miller_array = get_input_miller_array,
            align_datasets = align_datasets,
            extract_statistics = extract_dataset_statistics,
            copy_input_files = copy_input_files,
            write_output = write_initial_dataset_summary,
        )
        self.dataset_initialiser.map_scaling = config.params.maps.density_scaling
        self.dataset_initialiser.map_resolution_factor = config.params.maps.resolution_factor

    def process_objects(self, config):

        ###
        # Partitioners
        #

        from giant.mulch.partitioners import ResolutionShellTestTrainPartitioner
        self.partition_resolution_shells = ResolutionShellTestTrainPartitioner(
            shell_thickness = config.params.analysis.high_res_increment,
            high_resolution = config.params.analysis.high_res_upper_limit,
            low_resolution = config.params.analysis.high_res_lower_limit,
            test_train_partitioner = self.partition_test_train,
            min_train_datasets = config.params.statistical_maps.min_build_datasets,
            )

        ###
        # Map loading
        #

        from giant.mulch.transform.maps import (
            GetWarpedMapGrid,
            GetWarpedMapLoader,
            )

        self.get_map_grid = GetWarpedMapGrid(
            map_grid_spacing = float(config.params.maps.grid_spacing),
            outer_mask_radius = float(config.params.masks.outer_mask),
            inner_mask_radius = float(config.params.masks.inner_mask),
            symmetry_mask_radius = float(config.params.masks.inner_mask_symmetry),
            mask_pdb = config.params.masks.mask_pdb,
            align_mask_to_reference = bool(config.params.masks.align_mask_to_reference),
            create_grid_selection_string = config.params.masks.mask_selection_string,
            #mask_grid_selection_string = "not hetero",
            #partition_grid_selection_string = "pepnames and name CA and (altloc ' ' or altloc 'A')",
            processor = self.process_general,
            )

        self.get_warped_map_loader = GetWarpedMapLoader(
            get_map_grid = self.get_map_grid,
            processor = self.process_in_shell,
            )

        ###
        # Statistical model object
        #

        from pandda.analyse import (
            RunPanddaModel,
            GetPanddaModelProcessor,
            GetPanddaStatisticalModelFitter,
            GetPanddaStatisticalModelOutputter,
            GetPanddaDatasetEvaluator,
            GetPanddaDatasetEvaluatorOutputter,
            )

        from pandda.analyse.evaluate.maps import (
            GetStandardSparseMapWriter,
            )

        from pandda.analyse.events import (
            BasicPanddaFindEvents,
            BasicClusterFinder,
            ClusterFilterList,
            PeakAndSizeClusterFilter,
            GroupNearbyClustersFilter,
            ContactsClusterFilter,
            SymmetryClusterFilter,
            EventAnalyser,
            )

        # MASK NAMES
        load_map_mask_name = "outer"
        search_map_mask_name = "total"
        uncertainty_map_mask_name = "inner"

        get_pandda_model_task = GetPanddaModelProcessor(
            # Fit statistical model to input maps
            get_model_fitter = GetPanddaStatisticalModelFitter(
                processor = self.process_in_shell,
                fit_mu = config.params.statistical_maps.fit_mu,
                fit_sigma_adjusted = config.params.statistical_maps.fit_sigma_adjusted,
                get_map_grid = self.get_map_grid,
                load_map_mask_name = load_map_mask_name,
                sigma_uncertainty_map_mask_name = uncertainty_map_mask_name,
                ),
            # Find events, bdcs, etc
            get_dataset_evaluator = GetPanddaDatasetEvaluator(
                get_map_grid = self.get_map_grid,
                map_mask_name = load_map_mask_name,
                search_mask_name = search_map_mask_name,
                find_events = BasicPanddaFindEvents(
                    z_map_type = "adjusted+uncertainty",
                    find_clusters = BasicClusterFinder(
                        grid_clustering_cutoff = float(config.params.z_map_analysis.agglomerative_hierarchical.clustering_cutoff),
                        negative_values = bool(config.params.z_map_analysis.negative_values),
                        cluster_method = config.params.z_map_analysis.clustering_method,
                        contour_level = float(config.params.z_map_analysis.contour_level),
                        ),
                    filter_clusters = ClusterFilterList(
                        filters = [
                            PeakAndSizeClusterFilter(
                                min_blob_z_peak = float(
                                    config.params.z_map_analysis.min_blob_z_peak
                                    ),
                                min_blob_volume = float(
                                    config.params.z_map_analysis.min_blob_volume # angstroms
                                    ),
                                ),
                            ContactsClusterFilter(
                                distance_cutoff = 6.0, # angstroms
                                ),
                            GroupNearbyClustersFilter(
                                distance_cutoff = 5.0, # angstroms
                                ), 
                            SymmetryClusterFilter(
                                cluster_distance_cutoff = 8.0, # angstroms
                                contact_distance_cutoff = 6.0, # angstroms
                                ),
                            ],
                        ),
                    ),
                analyse_event = EventAnalyser(
                    max_bdc = config.params.background_correction.max_bdc,
                    min_bdc = config.params.background_correction.min_bdc,
                    increment = config.params.background_correction.increment,
                    output_multiplier = config.params.background_correction.output_multiplier,
                    output_directory = (
                        config.output.out_dir / "analysis" / "dataset_graphs" 
                        ),
                    ),
                processor = self.process_in_shell,
                ),
            # Writing functions
            get_model_outputter = GetPanddaStatisticalModelOutputter(
                output_dir = (
                    config.output.out_dir / "analysis"
                    ),
                processor = self.process_in_shell,
                ),
            get_dataset_outputter = GetPanddaDatasetEvaluatorOutputter(
                get_dataset_map_writer = GetStandardSparseMapWriter(
                    get_map_grid = self.get_map_grid,
                    mask_name = load_map_mask_name,
                    ),
                output_dir = (
                    config.output.out_dir / "analysis"
                    ),
                dataset_dir = (
                    config.output.out_dir / "processed_datasets"
                    ),
                processor = self.process_in_shell,
                output_requires_events = (config.output.output_maps_for == 'events'),
                ),
            )

        ###
        # Main processing function!
        #
        self.run_pandda_model = RunPanddaModel(
            get_pandda_model_task = get_pandda_model_task,
            partition_shells = self.partition_resolution_shells,
            output_dir = (
                config.output.out_dir / "shells"
                ),
            processor = self.process_shells,
            )
        #
        ###

        ###
        # Results function
        #
        from pandda.analyse.output import (
            PanddaResultsOutputter,
            MakePanddaResultsHtml,
            )
        #
        self.write_results = PanddaResultsOutputter(
            output_dir = (
                config.output.out_dir / "results"
                ),
            write_html = MakePanddaResultsHtml(
                output_directory = (
                    config.output.out_dir / "html"
                    ),
                ),
            )
        #
        ###

    def postprocess_objects(self, config):

        pass

        # AUTOBUILDING OBJECTS HERE


    def output_objects(self, config):

        from pandda.config import DumpConfigToJson
        self.dump_config_to_json = DumpConfigToJson(
            output_path = (
                config.output.out_dir / "params.json"
                ),
            )

        from pandda.output import DumpDictToJson
        self.dump_results_to_json = DumpDictToJson(
            output_path = (
                config.output.out_dir / "results.json"
                ),
            )

        from pandda.output.html import MakeMainPanddaHtmlPage
        self.make_html_output = MakeMainPanddaHtmlPage(
            output_directory = (
                config.output.out_dir / "html"
                ),
            dataset_html = str(self.dataset_initialiser.write_output.write_html.output_path),
            analyse_html = str(self.write_results.write_html.output_path),
            inspect_html = str(config.output.out_dir / "html" / "pandda_inspect.html"),
            )
