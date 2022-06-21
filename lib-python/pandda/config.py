import giant.logs as lg
logger = lg.getLogger(__name__)

import pathlib as pl

try:
    VERSION = pkg_resources.get_distribution("panddas").version
except:
    VERSION = '(developer -- see setup.py file)'


########################################################################################################################


class DumpConfigToJson(object):

    def __init__(self, output_path):

        self.output_path = pl.Path(output_path)

    def __call__(self, config):

        records = self.make_records(config)

        self.write_json(records=records)

        return self.output_path

    def make_records(self, config):

        records = {
            "data_dirs": str(pl.Path(config.input.data_dirs).absolute()),
            "out_dir": str(config.output.out_dir.absolute()),
        }

        return records

    def write_json(self, records):

        import json
        
        json_string = json.dumps(records, indent=2)

        with open(str(self.output_path), "w") as f:
            f.write(json_string)


########################################################################################################################


class InputError(Exception):
    pass


class ConfigurationError(Exception):
    pass


class Input(object):

    def __init__(self, config_obj):
        self.data_dirs = config_obj.data_dirs
        self.pdb_style = config_obj.pdb_style
        self.mtz_style = config_obj.mtz_style
        self.lig_style = config_obj.lig_style
        self.regex = InputRegex(config_obj.regex)
        self.filter = InputFilter(config_obj.filter)
        self.reference = InputReference(config_obj.reference)
        self.flags = InputFlags(config_obj.flags)


class InputRegex(object):

    def __init__(self, config_obj):
        self.dir_regex = config_obj.dir_regex
        self.pdb_regex = config_obj.pdb_regex
        self.mtz_regex = config_obj.mtz_regex
        # self.ligand_cif_regex = config_obj.ligand_cif_regex
        # self.ligand_pdb_regex = config_obj.ligand_pdb_regex


class InputFilter(object): 

    def __init__(self, config_obj):

        self.validate(config_obj)

        self.pdb = (
            pl.Path(config_obj.pdb) 
            if (config_obj.pdb is not None)
            else None
            )

    def validate(self, config_obj):

        if (config_obj.pdb is not None) and not pl.Path(config_obj.pdb).exists():
            raise InputError(
                "filter pdb does not exist: {}".format(
                    config_obj.pdb
                    )
                )


class InputReference(object):

    def __init__(self, config_obj):

        self.validate(config_obj)

        self.pdb = config_obj.pdb
        self.mtz = config_obj.mtz
        #self.structure_factors = config_obj.structure_factors.split(',')

    def validate(self, config_obj):

        if (config_obj.pdb is not None) and not pl.Path(config_obj.pdb).exists():
            raise InputError(
                "reference pdb does not exist: {}".format(
                    config_obj.pdb
                    )
                )

        if (config_obj.mtz is not None) and not pl.Path(config_obj.mtz).exists():
            raise InputError(
                "reference mtz does not exist: {}".format(
                    config_obj.mtz
                    )
                )


class InputFlags(object):

    def __init__(self, config_obj):

        self.validate(config_obj)

        self.ignore_datasets = (
            config_obj.ignore_datasets.split(',')
            if config_obj.ignore_datasets
            else None
            )

        self.only_datasets = (
            config_obj.only_datasets.split(',')
            if config_obj.only_datasets
            else None
            )

        self.test = (
            config_obj.test.split(',')
            if config_obj.test 
            else None
            )

        self.train = (
            config_obj.train.split(',')
            if config_obj.train
            else config_obj.ground_state_datasets.split(',')
            if config_obj.ground_state_datasets
            else None
            )

        self.not_test = (
            config_obj.not_test.split(',')
            if config_obj.not_test
            else config_obj.exclude_from_z_map_analysis.split(',')
            if config_obj.exclude_from_z_map_analysis
            else None
            )
        
        self.not_train = (
            config_obj.not_train.split(',')
            if config_obj.not_train
            else config_obj.exclude_from_characterisation.split(',')
            if config_obj.exclude_from_characterisation
            else None
            )

    def validate(self, config_obj):

        if (config_obj.ground_state_datasets and config_obj.train):
            raise ValueError(
                "Can supply ground_state_datasets OR train, but not both."
                )

        if (config_obj.exclude_from_z_map_analysis and config_obj.not_test):
            raise ValueError(
                "Can supply exclude_from_z_map_analysis OR not_test, but not both."
                )

        if (config_obj.exclude_from_characterisation and config_obj.not_train):
            raise ValueError(
                "Can supply exclude_from_characterisation OR not_train, but not both."
                )


########################################################################################################################


class Output(object):

    def __init__(self, config_obj):

        self.validate(config_obj)

        self.out_dir = pl.Path(config_obj.out_dir)
        self.overwrite = config_obj.overwrite
        self.dataset_prefix = config_obj.dataset_prefix
        self.output_maps_for = config_obj.output_maps_for

    def validate(self, config_obj):

        if pl.Path(config_obj.out_dir).exists():
            if not pl.Path(config_obj.out_dir).is_dir():
                raise InputError(
                    "Output directory already exists but is not a directory: {}".format(
                        config_obj.out_dir
                        )
                    )
            if config_obj.overwrite: 
                import shutil
                shutil.rmtree(config_obj.out_dir)
            else:
                raise InputError(
                    "Output directory already exists: {}".format(
                        config_obj.out_dir,
                        )
                    )
        else: 
            pl.Path(config_obj.out_dir).mkdir()


########################################################################################################################


class Params(object):

    def __init__(self, config_obj):
        self.analysis = Analysis(config_obj.analysis)
        self.diffraction_data = DiffractionData(config_obj.diffraction_data)
        self.alignment = Alignment(config_obj.alignment)
        self.filtering = Filtering(config_obj.filtering)
        self.excluding = Excluding(config_obj.excluding)
        self.maps = Maps(config_obj.maps)
        self.masks = Masks(config_obj.masks)
        self.statistical_maps = StatisticalMaps(config_obj.statistical_maps)
        self.z_map_analysis = ZMapAnalysis(config_obj.z_map_analysis)
        self.background_correction = BackgroundCorrection(config_obj.background_correction)


class Analysis(object):

    def __init__(self, config_obj):

        self.high_res_upper_limit = float(config_obj.high_res_upper_limit)
        self.high_res_lower_limit = float(config_obj.high_res_lower_limit)
        self.high_res_increment = float(config_obj.high_res_increment)


class DiffractionData(object):

    def __init__(self, config_obj):

        self.validate(config_obj)

        self.structure_factors = [
            tuple(sf.split(',')) 
            for sf in config_obj.structure_factors
            ]
        self.checks = DiffractionDataChecks(config_obj.checks)
        self.scaling = DiffractionDataScaling(config_obj.scaling)

    def validate(self, config_obj):

        for sf in config_obj.structure_factors:
            if (',' not in sf):
                raise InputError("structure_factors must be F,PHI pairs separated by a comma.")


class DiffractionDataChecks(object):

    def __init__(self, config_obj):

        self.validate(config_obj)

        self.low_resolution_completeness = config_obj.low_resolution_completeness
        self.all_data_are_valid_values = config_obj.all_data_are_valid_values

    def validate(self, config_obj):

        if config_obj.low_resolution_completeness <= 0.0: 
            raise InputError("low_resolution_completeness must be greater than 0.0")


class DiffractionDataScaling(object):

    def __init__(self, config_obj):

        self.apply_b_factor_scaling = config_obj.apply_b_factor_scaling


class Alignment(object):

    def __init__(self, config_obj):
        self.method = config_obj.method


class Filtering(object):

    def __init__(self, config_obj):
        self.max_rfree = config_obj.max_rfree
        self.same_space_group_only = config_obj.same_space_group_only
        self.similar_models_only = config_obj.similar_models_only


class Excluding(object):

    def __init__(self, config_obj):
        self.max_wilson_plot_z_score = config_obj.max_wilson_plot_z_score


class Maps(object):

    def __init__(self, config_obj):
        self.resolution_factor = config_obj.resolution_factor
        self.grid_spacing = config_obj.grid_spacing
        self.density_scaling = config_obj.density_scaling


class Masks(object):

    def __init__(self, config_obj):
        self.mask_pdb = config_obj.mask_pdb
        self.mask_selection_string = config_obj.mask_selection_string
        self.align_mask_to_reference = config_obj.align_mask_to_reference
        self.outer_mask = config_obj.outer_mask
        self.inner_mask = config_obj.inner_mask
        self.inner_mask_symmetry = config_obj.inner_mask_symmetry


class ZMapAnalysis(object):

    def __init__(self, config_obj):
        self.clustering_method = config_obj.clustering_method
        self.contour_level = config_obj.contour_level
        self.negative_values = config_obj.negative_values

        self.min_blob_volume = config_obj.min_blob_volume
        self.min_blob_z_peak = config_obj.min_blob_z_peak
        self.masks = ZMasks(config_obj.masks)
        self.agglomerative_hierarchical = AgglomerativeHierarchical(config_obj.agglomerative_hierarchical)


class ZMasks(object):

    def __init__(self, config_obj):
        self.selection_string = config_obj.selection_string
        self.outer_mask = config_obj.outer_mask
        self.inner_mask = config_obj.inner_mask


class BackgroundCorrection(object):

    def __init__(self, config_obj):
        self.max_bdc = config_obj.max_bdc
        self.min_bdc = config_obj.min_bdc
        self.increment = config_obj.increment
        self.output_multiplier = config_obj.output_multiplier


class AgglomerativeHierarchical(object):

    def __init__(self, config_obj):
        self.clustering_cutoff = config_obj.clustering_cutoff


class StatisticalMaps(object):

    def __init__(self, config_obj):
        self.min_build_datasets = config_obj.min_build_datasets
        self.max_build_datasets = config_obj.max_build_datasets
        self.fit_mu = config_obj.fit_mu
        self.fit_sigma_adjusted = config_obj.fit_sigma_adjusted


########################################################################################################################


class Results(object):

    def __init__(self, config_obj):
        self.events = Events(config_obj.events)


class Events(object):

    def __init__(self, config_obj):
        self.order_by = config_obj.order_by


########################################################################################################################


class Processing(object):

    def __init__(self, config_obj):
        
        self.cpus = config_obj.cpus
        self.backend = config_obj.backend
        
        self.process_shells = config_obj.process_shells
        self.remote_nodes = config_obj.remote_nodes
        self.cpus_per_remote_node = config_obj.cpus_per_remote_node
        
        self.h_vmem = config_obj.h_vmem
        self.m_mem_free = config_obj.m_mem_free


########################################################################################################################


class Settings(object):

    def __init__(self, config_obj):
        
        self.verbose = config_obj.verbose
        self.plotting_backend = config_obj.plotting.backend


########################################################################################################################


class Autobuilding(object):

    def __init__(self, config_obj):
        self.autobuild = config_obj.autobuild


########################################################################################################################


class Config(object):

    def __init__(self, config_obj):
        self.settings = Settings(config_obj.settings)

        self.input = Input(config_obj.pandda.input)
        self.output = Output(config_obj.pandda.output)
        self.params = Params(config_obj.pandda.params)
        self.processing = Processing(config_obj.pandda.processing)
        self.results = Results(config_obj.pandda.results)
        #self.autobuilding = Autobuilding(config_obj.pandda.autobuilding)
