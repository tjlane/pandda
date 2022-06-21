import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy as np

from giant.mulch.partitioners import (
    DatasetKeyPartitioner,
    )

from pandda.utils import (
    merge_dicts,
    )


class DatasetInitialiser(object):

    map_scaling = "rmsd"
    map_resolution_factor = 0.25

    def __init__(self,
        dataloader,
        partitioner,
        initial_filter,
        select_reference,
        filter_with_reference,
        check_input_data,
        calculate_scalings,
        get_input_miller_array,
        align_datasets,
        extract_statistics,
        copy_input_files, # change to "copy input files"
        write_output,
        ):

        self.dataloader = dataloader
        self.initial_filter = initial_filter
        self.select_reference = select_reference
        self.filter_with_reference = filter_with_reference
        self.check_input_data = check_input_data
        self.calculate_scalings = calculate_scalings
        self.get_input_miller_array = get_input_miller_array
        self.align_datasets = align_datasets
        self.extract_statistics = extract_statistics
        self.copy_input_files = copy_input_files
        self.write_output = write_output

        # Record list of filters
        self.filters = []

        # Total partitioner from the class
        self.partitioner = partitioner

    def __call__(self):

        ###
        #
        # Load and filter datasets
        #
        ###

        logger.heading('Loading datasets')
        mcd = self.dataloader()
        logger(str(self.dataloader))

        if (self.initial_filter is not None):
            logger.subheading('Filtering input datasets')
            mcd = self.initial_filter(
                mcd = mcd,
                )
            self.filters.append(self.initial_filter)
            logger(str(self.initial_filter))

        logger.subheading('Selecting refererence dataset')
        reference_dataset = self.select_reference(
            mcd = mcd,
            )
        mcd.set_reference_dataset(reference_dataset)

        if (self.filter_with_reference is not None):
            logger.subheading('Filtering by reference structure')
            mcd = self.filter_with_reference(
                mcd = mcd,
                reference_dataset = reference_dataset,
                )
            self.filters.append(self.filter_with_reference)
            logger(str(self.filter_with_reference))

        if (self.check_input_data is not None):
            logger.subheading('Checking input datasets')
            self.check_input_data(
                mcd = mcd,
                )

        if (self.calculate_scalings is not None):
            logger.subheading('Scaling input datasets to reference')
            scaling_objects = self.calculate_scalings(
                mcd = mcd,
                reference_dataset = reference_dataset,
                )
            scaling_filter =self.calculate_scalings.as_filter()
            mcd = scaling_filter(mcd)
            self.filters.append(scaling_filter)
            logger(str(scaling_filter))
        else:
            scaling_objects = {}

        if (self.align_datasets is not None):
            logger.subheading('Aligning input structures to reference')
            alignment_objects = self.align_datasets(
                mcd = mcd, 
                reference_dataset = reference_dataset,
                )
            alignments_filter = self.align_datasets.as_filter()
            mcd = alignments_filter(mcd)
            self.filters.append(alignments_filter)
            logger(str(alignments_filter))
        else:
            alignment_objects = {}

        ###
        #
        # Calculate output statistics and run output functions
        #
        ###

        if (self.extract_statistics is not None):
            logger.subheading('Extracting data statistics')
            dataset_statistics = self.extract_statistics(
                mcd = mcd,
                get_miller_array = self.get_input_miller_array,
                get_scaling_object = lambda d: scaling_objects.get(d.tag, None),
                )
            logger(str(dataset_statistics))
            # Add partitions from this object to global partitions
            self.partitioner.update(
                partitions_dict = dataset_statistics.classifications,
                )
            logger(str(self.partitioner))

        ###
        #
        # Construct output objects
        #
        ###

        from giant.mulch.xray import DataGetterPrepper, MapGetterPrepper
        self.data_getters = {
            dtag: DataGetterPrepper(
                get_data = self.get_input_miller_array,
                truncate_data = None, # added later
                scale_data = scaling_objects.get(dtag, None),
                )
            for dtag, dataset in mcd.datasets.items()
            }

        self.map_getters = {
            dtag: MapGetterPrepper(
                get_miller_array = self.data_getters[dtag],
                map_scaling = self.map_scaling, 
                map_resolution_factor = self.map_resolution_factor,
                )
            for dtag, dataset in mcd.datasets.items()
            }

        self.reference_map_getter = MapGetterPrepper(
            get_miller_array = DataGetterPrepper(
                get_data = self.get_input_miller_array,
                truncate_data = None, # added later if ever
                scale_data = None, # no scaling for the reference dataset
                ),
            map_scaling = self.map_scaling, 
            map_resolution_factor = self.map_resolution_factor,
            )

        ###
        #
        # Write output
        #
        ###

        output_files = {
            'reference_files' : {
                'reference_model' : reference_dataset.model.filename,
                'reference_data' : reference_dataset.data.filename,
                },
        }

        if (self.copy_input_files is not None):
            of = self.copy_input_files(
                mcd = mcd,
                )
            merge_dicts(
                master_dict = output_files,
                merge_dict = {'dataset_files' : of},
                )

        if (self.write_output is not None):
            logger.subheading('Writing initial dataset files')
            of = self.write_output(
                mcd = mcd,
                dataset_statistics = dataset_statistics,
                dataloader = self.dataloader,
                filters = self.filters,
                partitions = self.partitioner.as_dict(),
                data_getters = self.data_getters,
                reference_dataset = reference_dataset,
                get_reference_data = self.get_input_miller_array,
                )
            merge_dicts(
                master_dict = output_files,
                merge_dict = of,
                )

        logger.subheading('Input dataset summary')
        logger(str(mcd))

        self.output_files = output_files
        
        return mcd

    def __str__(self):
        return ""

