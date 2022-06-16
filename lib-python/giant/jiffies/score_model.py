import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys, copy, itertools

import numpy as np
import pathlib as pl

from giant.mulch.labelling import (
    PathLabeller,
    )

from giant.exceptions import (
    Sorry, Failure,
    )

from giant.processors import (
    ProcessorJoblib,
    )

from giant.phil import (
    log_running_parameters,
    )

from giant.plot import (
    setup_plotting,
    )

from giant.mulch.dataset import (
    CrystallographicDataset,
    )

from giant.validation.score_residues import (
    ScoreModelSingle,
    ScoreModelMultiple,
    GetInterestingResnames,
    WriteScores,
    ValidationRadarPlot
    )

setup_plotting()

#################################

PROGRAM = 'giant.score_model'
DESCRIPTION = """
    A tool to quickly score residues (e.g. ligands) against crystallographic electron density.

    Single Structures:

    1) Simple usage (for ligand called LIG or UNL):
        > giant.score_model <filename>.pdb <filename>.mtz

    2) Score particular residues in the file (residues XX1, XX2 and XX3)
        > giant.score_model ... resnames=XX1,XX2,XX3

    3) Define a "reference" model and data to compare the input model to. Providing the MTZ is optional
        > giant.score_model ... ref_pdb=<filename>.pdb ref_mtz=<filename>.mtz
        > giant.score_model ... ref_pdb=<filename>.pdb

    Multiple Structures:

    1) Simple usage:
        Mtz files can be provided explicitly:
        > giant.score_model filename1.pdb filename1.mtz filename2.pdb filename2.mtz filename3.pdb filename3.mtz
        or can be omitted when the same same as the pdb files (e.g. filename1.mtz, filename2.mtz, etc)
        > giant.score_model filename1.pdb filename2.pdb filename3.pdb

    See giant.score_model for more information.
"""

#######################################

blank_arg_prepend = {
    '.pdb' : 'pdb=',
    '.mtz' : 'mtz=',
     None  : 'directory=',
}

import libtbx.phil
master_phil = libtbx.phil.parse("""
input {
    from_files {
        pdb = None
            .type = str
            .multiple = True
        mtz = None
            .type = str
            .multiple = True
        reference_pdb = None
            .type = str
            .multiple = True
        reference_mtz = None
            .type = str
            .multiple = True
    }
    from_directories {
        directory = None
            .type = str
            .multiple = True
        pdb_style = None
            .type = str
        mtz_style = None
            .type = str
        reference_pdb_style = None
            .type = str
        reference_mtz_style = None
            .type = str
    }
    label = None
        .help = "Provide labels for the input structure(s)"
        .type = str
        .multiple = False
    labelling = *automatic basename filename foldername
        .help = "If labels are not provided, labels will be derived from the filename of the input structure"
        .type = choice(multi=False)
}
output {
    out_dir = score_model
        .type = path
    log = None
        .type = path
}
options {
    ignore_common_molecules = True
        .type = bool
    include_resname = None
        .type = str
        .multiple = True
    ignore_resname = None
        .type = str
        .multiple = True
    edstats_f_label = None
        .help = 'Column label for experimental amplitudes in input mtz files. If left blank will try to guess the labels.'
        .type = str
        .multiple = False
}
plot {
    remove_blank_entries = False
        .type = bool
    parameters {
        rscc {
            title = '\n\nModel\nQuality\n(RSCC)'
                .type = str
            axis_min = 0.60
                .type = float
            axis_max = 0.85
                .type = float
            axis_invert = True
                .type = bool
        }
        rszd {
            title = 'Model\nAccuracy\n(RSZD)'
                .type = str
            axis_min = 1.50
                .type = float
            axis_max = 4.00
                .type = float
            axis_invert = False
                .type = bool
        }
        rszo {
            title = 'Density\nPrecision\n(RSZO/OCC)'
                .type = str
            axis_min = 0.00
                .type = float
            axis_max = 2.00
                .type = float
            axis_invert = True
                .type = bool
        }
        b_factor_ratio {
            title = 'B-Factor\nStability\n(B-factor Ratio)'
                .type = str
            axis_min = 1.00
                .type = float
            axis_max = 3.00
                .type = float
            axis_invert = False
                .type = bool
        }
        rmsd {
            title = 'Coordinate\nStability\n(RMSD)'
                .type = str
            axis_min = 0.00
                .type = float
            axis_max = 1.50
                .type = float
            axis_invert = False
                .type = bool
        }
    }
}
include scope giant.phil.settings_phil
""",
process_includes=True,
)


#######################################


class validate_params(object):

    def __init__(self, params):

        self.validate_input(params)

    def validate_input(self, params):

        self.validate_structures(params)
        self.validate_reference_files(params)
        self.validate_labels(params)

    def validate_structures(self, params):

        if (params.input.from_files.pdb and params.input.from_directories.directory):
            raise Sorry('Must provide input structures OR input directories')

        if (not params.input.from_files.pdb) and (not params.input.from_directories.directory):
            raise Sorry('No input structures provided (provide structures OR input directories')

        if (params.input.from_files.pdb):

            self.validate_input_files(
                params.input.from_files
                )

        else:

            self.validate_input_from_directories(
                params.input.from_directories
                )

    def validate_input_files(self, f_params):

        for p in f_params.pdb:

            assert pl.Path(p).exists()

        if (not f_params.mtz):

            logger(
                'No MTZ files provided. Looking for files with the same naming as pdb files: \n\t{}'.format(
                    '\n\t'.join([
                        '{pdb} -> {mtz}'.format(
                            pdb = p,
                            mtz = str(pl.Path(p).with_suffix('.mtz')),
                            )
                        for p in f_params.pdb
                        ]),
                )
            )

            mtz_files = []

            for p in f_params.pdb:

                m = pl.Path(p).with_suffix('.mtz')

                if not m.exists():
                    raise ValueError('No MTZs provided and MTZ not found: {}'.format(str(m)))

                f_params.mtz.append(str(m))

        else:

            for m in f_params.mtz:

                assert pl.Path(m).exists()

        assert len(f_params.pdb) == len(f_params.mtz)

    def validate_input_from_directories(self, f_params):

        assert f_params.directory
        assert f_params.pdb_style is not None

        if f_params.mtz_style is None:

            f_params.mtz_style = (
                f_params.pdb_style.rsplit('.pdb')[0] + '.mtz'
                )

            logger('Setting mtz_style is automatic value: {}'.format(f_params.mtz_style))

        dir_paths = list(map(pl.Path, f_params.directory))

        for d in dir_paths:

            pdbs = list(d.glob(f_params.pdb_style))
            mtzs = list(d.glob(f_params.mtz_style))

            if len(pdbs) == 0:
                raise ValueError(
                    'No PDB files found in directory for pdb_style: {}'.format(
                        str(d),
                        )
                    )

            if len(pdbs) > 1:
                raise ValueError(
                    'More than one PDB file in directory {} matches pdb_style: \n\t{}'.format(
                        str(d),
                        str(pdbs),
                        )
                    )

            if len(mtzs) == 0:
                raise ValueError(
                    'No MTZ files found in directory for mtz_style: {}'.format(
                        str(d)
                        )
                    )
            if len(mtzs) > 1:
                raise ValueError(
                    'More than one MTZ file in directory {} matches mtz_style: \n\t{}'.format(
                        str(d),
                        str(mtzs),
                        )
                    )

    def validate_reference_files(self, params):

        if (params.input.from_files.reference_pdb):

            for p in params.input.from_files.reference_pdb:

                assert pl.Path(p).exists()

        if (params.input.from_files.reference_mtz):

            assert (
                len(params.input.from_files.reference_pdb) ==
                len(params.input.from_files.reference_mtz)
                )

            for m in params.input.from_files.reference_mtz:

                assert pl.Path(m).exists()

        elif (params.input.from_files.reference_pdb):

            for p in params.input.from_files.reference_pdb:

                m = pl.Path(p).with_suffix('.mtz')

                if m.exists():

                    params.input.from_files.reference_mtz.append(str(m))

                else:

                    params.input.from_files.reference_mtz.append(None)

    def validate_labels(self, params):

        if (params.input.label) and (params.input.pdb):

            assert len(params.input.label) == len(params.input.pdb)

        elif (params.input.label) and (params.input.directory):

            assert len(params.input.label) == len(params.input.directory)


class Config(object):

    def __init__(self, params):

        self.params = params

    def setup_dir(self):

        p = self.params

        out_dir_path = pl.Path(p.output.out_dir)
        if not out_dir_path.exists():
            out_dir_path.mkdir(parents=True)

        return out_dir_path

    def set_get_log(self):

        p = self.params

        if p.output.log:
            pass
        else:
            p.output.log = str(
                pl.Path(p.output.out_dir) / 'score_model.log'
                )

        return p.output.log

    def get_datasets(self):

        logger.subheading('Loading Datasets')

        p = self.params

        # Apply labels

        labels = (
            list(p.input.label)
            if p.input.label is not None
            else None
            )

        #

        datasets = []

        if p.input.from_files.pdb:

            if labels:

                assert len(labels) == len(p.input.from_files.pdb)

            else:

                labelling = (
                    PathLabeller('basename')
                    if
                    p.input.labelling == 'automatic'
                    else
                    PathLabeller(p.input.labelling)
                    )

            for pdb_path, mtz_path in zip(
                p.input.from_files.pdb,
                p.input.from_files.mtz,
                ):

                d = CrystallographicDataset.from_file(
                    model_filename = pdb_path,
                    data_filename = mtz_path,
                    )

                if labels:
                    d.label(tag=labels.pop(0))
                else:
                    d.label(tag=labelling(pdb_path))

                datasets.append(d)

        elif p.input.from_directories.directory:

            if labels:

                assert len(labels) == len(p.input.from_directories.directory)

            else:

                labelling = (
                    PathLabeller('foldername')
                    if
                    p.input.labelling == 'automatic'
                    else
                    PathLabeller(p.input.labelling)
                    )

            for dir_path in map(pl.Path, p.input.from_directories.directory):

                pdb_path = next(dir_path.glob(p.input.from_directories.pdb_style))
                mtz_path = next(dir_path.glob(p.input.from_directories.mtz_style))

                d = CrystallographicDataset.from_file(
                    model_filename = str(pdb_path),
                    data_filename = str(mtz_path),
                    )

                if labels:
                    d.label(tag=labels.pop(0))
                else:
                    d.label(tag=labelling(pdb_path))

                datasets.append(d)

        else:

            raise NotImplementedError()

        ##

        logger('\n> Loaded datasets:\n')

        for d in datasets:

            logger(str(d))

        logger('\n\n> {} datasets loaded'.format(len(datasets)))

        ##

        return datasets

    def get_reference_datasets(self):

        logger.subheading('Loading Reference Datasets')

        p = self.params

        datasets = []

        if p.input.from_files.reference_pdb:

            for pdb_path, mtz_path in zip(
                p.input.from_files.reference_pdb,
                p.input.from_files.reference_mtz,
                ):

                d = CrystallographicDataset.from_file(
                    model_filename = pdb_path,
                    data_filename = mtz_path,
                    )

                datasets.append(d)

        elif p.input.from_directories.reference_pdb_style:

            for dir_path in map(pl.Path, p.input.from_directories.directory):

                try:
                    pdb_path = str(
                        next(dir_path.glob(p.input.from_directories.reference_pdb_style))
                        )
                except:
                    logger(
                        'No reference PDB found in {}'.format(
                            str(dir_path)
                            )
                        )
                    pdb_path = None

                if pdb_path is None:
                    datasets.append(None)
                    continue

                try:
                    mtz_path = str(
                        next(dir_path.glob(p.input.from_directories.reference_mtz_style))
                        )
                except:
                    logger(
                        'No reference MTZ found for {}'.format(
                            str(dir_path)
                            )
                        )
                    mtz_path = None

                d = CrystallographicDataset.from_file(
                    model_filename = pdb_path,
                    data_filename = mtz_path,
                    )

                datasets.append(d)

        else:

            logger('No reference datasets loaded')

            return None

        ##

        logger('\n> Loaded reference datasets:\n')

        for d in datasets:

            logger(str(d))

        logger('\n\n> {} reference datasets loaded'.format(len(datasets)))

        ##

        return datasets

    def get_plot_parameters(self):

        p = self.params.plot.parameters

        d = {
            'axis_params' : {
                'rscc' : {
                    'title' : p.rscc.title,
                    'axis_min' : p.rscc.axis_min,
                    'axis_max' : p.rscc.axis_max,
                    'axis_invert' : p.rscc.axis_invert,
                },
                'rszd' : {
                    'title' : p.rszd.title,
                    'axis_min' : p.rszd.axis_min,
                    'axis_max' : p.rszd.axis_max,
                    'axis_invert' : p.rszd.axis_invert,
                },
                'rszo' : {
                    'title' : p.rszo.title,
                    'axis_min' : p.rszo.axis_min,
                    'axis_max' : p.rszo.axis_max,
                    'axis_invert' : p.rszo.axis_invert,
                },
                'b_factor_ratio' : {
                    'title' : p.b_factor_ratio.title,
                    'axis_min' : p.b_factor_ratio.axis_min,
                    'axis_max' : p.b_factor_ratio.axis_max,
                    'axis_invert' : p.b_factor_ratio.axis_invert,
                },
                'rmsd' : {
                    'title' : p.rmsd.title,
                    'axis_min' : p.rmsd.axis_min,
                    'axis_max' : p.rmsd.axis_max,
                    'axis_invert' : p.rmsd.axis_invert,
                },
            },
        }

        return d


def run(params):

    config = Config(params)

    config.setup_dir()

    logger = lg.setup_logging(
        name = __name__,
        log_file = config.set_get_log(),
        debug = bool(params.settings.verbose),
        )

    logger.heading(
        'Validating input parameters and input files'
        )

    validate_params(params)

    log_running_parameters(
        params = params,
        master_phil = master_phil,
        logger = logger,
    )

    logger.heading(
        'Setting up...'
        )

    score_models = ScoreModelMultiple(
        score_model = ScoreModelSingle(
            get_interesting_resnames = GetInterestingResnames(
                ignore_common_molecules = params.options.ignore_common_molecules,
                include_resnames_list = params.options.include_resname,
                ignore_resnames_list = params.options.ignore_resname,
                ),
            ),
        processor = ProcessorJoblib(
            n_cpus = params.settings.cpus,
            ),
        )

    write_output = WriteScores(
        make_radar_plot = ValidationRadarPlot(
            plot_defaults = config.get_plot_parameters(),
            ),
        )

    logger.heading('Loading Datasets')

    datasets = config.get_datasets()
    reference_datasets = config.get_reference_datasets()

    logger.heading('Scoring Models')

    results_df = score_models(
        dataset_list = datasets,
        reference_datasets = reference_datasets,
        )

    logger.heading('Writing output')

    output_dict = write_output(
        residue_scores_df = results_df,
        dir_path = pl.Path(params.output.out_dir),
        )

    if False:

        logger.subheading('Collating all images and outputting to output HTML...')

        # Get template to be filled in
        from giant.html import HTML_ENV
        template = HTML_ENV.get_template('summary_page.html')

        # Output directory (for relative symlinks)
        out_dir = os.path.abspath(os.path.dirname(summary_file))

        # ===========================================================>
        # Construct the data object to populate the template
        output_data = {}
        output_data['header'] = 'Residue Score Summaries'
        output_data['title'] = 'Residue Score Summaries'
        output_data['introduction'] = 'Model Quality and Validation checks.'
        # ===========================================================>
        # Header Images
        output_data['small_images'] = []
        for img in (list(combined_images.values()) + list(separate_images.values())):
            output_data['small_images'].append({
                'path': './'+os.path.relpath(path=img, start=out_dir),
                'title': 'Scores for {}'.format(os.path.splitext(os.path.basename(img))[0]),
            })
        # ===========================================================>
        # Write Output
        with open(summary_file, 'w') as out_html:
            out_html.write(template.render(output_data))

    logger.heading('finished normally!')

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
        program             = PROGRAM,
        description         = DESCRIPTION,
    )
