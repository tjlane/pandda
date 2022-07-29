import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys, tempfile, collections
import pandas as pd

from giant.exceptions import Sorry, Failure

_EdstatsGlobalScores = collections.namedtuple(
    typename = 'EdstatsGlobalScores',
    field_names = [
        'mean_biso',
        'rms_z_scores',
    ],
)

class EdstatsGlobalScores(_EdstatsGlobalScores):

    def __str__(self):

        s_ = (
            "EdstatsGlobalScores: \n"
            "\tMean B (iso): {mean_biso}\n"
            "\tRMS Z-scores: \n"
            "\t\t{rms_z_scores}\n"
        ).format(
            mean_biso = self.mean_biso,
            rms_z_scores = '\n\t\t'.join([
                '{k} : {v}'.format(
                    k = k,
                    v = v,
                    )
                for k,v in sorted(
                    self.rms_z_scores.items()
                    )
                ]),
            )

        return s_

_EdstatsResidueScores = collections.namedtuple(
    typename = 'EdstatsResidueScores',
    field_names = [
        'rscc',
        'rsr',
        'rszo',
        'rszd',
        'average_b',
    ],
)

class EdstatsResidueScores(_EdstatsResidueScores):

    def __str__(self):

        s_ = (
            "EdstatsResidueScores: \n"
            "\tMean B (iso): {mean_biso}\n"
            "\tRMS Z-scores: \n"
            "\t\t{rms_z_scores}\n"
        ).format(
            rscc = self.rscc,
            rsr = self.rsr,
            rszo = self.rszo,
            rszd = self.rszd,
            average_b = self.average_b,
            )

        return s_


class EdstatsResults(object):

    def __init__(self,
        results_table,
        global_scores = None,
        ):

        self.results_table = (
            results_table
            )

        self.global_scores = (
            global_scores
            )

    def __str__(self):

        s_ = (
            "EdstatsResults: \n"
            "\tGlobal Scores: \n\t\t{global_scores}\n"
            "\tResults Table: \n\t\t{n_results} residues scored\n"
            ).format(
            global_scores = (
                str(self.global_scores).replace('\n','\n\t\t')
                ),
            n_results = (
                len(self.results_table)
                )
            )

        return s_

    def get_edstats_label(self, residue_group):

        rg = residue_group

        model_no = rg.parent().parent().id.strip()
        chain_id = rg.parent().id.strip()
        res_id = rg.resid().strip()

        if model_no == '':
            model_no = '1'

        if chain_id == '':
            chain_id = '_'

        return (model_no, chain_id, res_id)

    def get(self, residue_group=None, residue_group_label=None):

        if [residue_group, residue_group_label].count(None) != 1:
            raise ValueError('Provide EITHER residue_group OR residue_group_label')

        # Shortcuts
        rg = residue_group
        rg_label = residue_group_label

        if (rg_label is None):
            rg_label = self.get_edstats_label(residue_group=rg)

        # Check label associated with a record
        if not self.results_table.index.contains(rg_label):
            raise KeyError('The given label does not have a corresponding result: {}'.format(rg_label))

        # Extract residue scores
        ed_scores = self.results_table.loc[rg_label]

        return EdstatsResidueScores(
            rscc = ed_scores['CCSa'],
            rsr = ed_scores['Ra'],
            rszo = ed_scores['ZOa'],
            rszd = ed_scores['ZDa'],
            average_b = ed_scores['BAa'],
            )


class ParseEdstatsLog(object):

    RMS_DEFAULTS = {
        'atm_all_ZD+' : None,
        'atm_all_ZD-' : None,
        'atm_all_ZO' : None,
        'atm_main_ZD+' : None,
        'atm_main_ZD-' : None,
        'atm_main_ZO' : None,
        'atm_side_ZD+' : None,
        'atm_side_ZD-' : None,
        'atm_side_ZO' : None,
        'het_all_ZD+' : None,
        'het_all_ZD-' : None,
        'het_all_ZO' : None,
        'het_main_ZD+' : None,
        'het_main_ZD-' : None,
        'het_main_ZO' : None,
        'het_side_ZD+' : None,
        'het_side_ZD-' : None,
        'het_side_ZO' : None,
    }

    def __init__(self):

        pass

    def __call__(self, log_text):

        output = [line for line in log_text.split('\n') if line]

        rms_scores = dict(self.RMS_DEFAULTS)
        mean_biso = None

        for i, line in enumerate(output):

            if line.startswith('Overall mean Biso:'):

                mean_biso = float(line.replace('Overall mean Biso:',''))

            elif line.startswith('RMS Z-scores for protein residues:'):

                assert ['Main','Side','All'] == output[i+1].strip().split(), (
                    'COLUMN HEADINGS ARE NOT WHAT IS EXPECTED! {!s}'.format(line)
                    )

                for incr in [2,3,4]:

                    score, main_val, side_val, all_val = output[i+incr].strip().split()

                    rms_scores['_'.join(['atm','main',score])] = float(main_val)
                    rms_scores['_'.join(['atm','side',score])] = float(side_val)
                    rms_scores['_'.join(['atm','all',score])] = float(all_val)

            elif line.startswith('RMS Z-scores for hetero residues:'):

                assert ['Main','Side','All'] == output[i+1].strip().split(), (
                    'COLUMN HEADINGS ARE NOT WHAT IS EXPECTED! {!s}'.format(line)
                    )

                for incr in [2,3,4]:

                    score, main_val, side_val, all_val = output[i+incr].strip().split()

                    rms_scores['_'.join(['het','main',score])] = float(main_val)
                    rms_scores['_'.join(['het','side',score])] = float(side_val)
                    rms_scores['_'.join(['het','all',score])] = float(all_val)

        return EdstatsGlobalScores(
            mean_biso = mean_biso,
            rms_z_scores = rms_scores,
            )


class CalculateEdstatsScores(object):

    _parse_log = ParseEdstatsLog()

    # Output columns in order
    EDSTATS_COLUMNS = [
        'RT', 'CI', 'RN',
        'BAm', 'NPm', 'Rm', 'RGm', 'SRGm', 'CCSm', 'CCPm', 'ZCCPm', 'ZOm', 'ZDm', 'ZD-m', 'ZD+m',
        'BAs', 'NPs', 'Rs', 'RGs', 'SRGs', 'CCSs', 'CCPs', 'ZCCPs', 'ZOs', 'ZDs', 'ZD-s', 'ZD+s',
        'BAa', 'NPa', 'Ra', 'RGa', 'SRGa', 'CCSa', 'CCPa', 'ZCCPa', 'ZOa', 'ZDa', 'ZD-a', 'ZD+a',
        'MN', 'CP', 'NR'
    ]

    # Output datatypes
    EDSTATS_DTYPES = {
        'RT' : str, # Keep as string for use as labels
        'CI' : str, # Keep as string for use as labels
        'RN' : str, # Keep as string for use as labels
        'MN' : str, # Keep as string for use as labels
        'CP' : str, # Keep as string for use as labels
        'NR' : str, # Keep as string for use as labels
        'BAm' : float, 'NPm' : pd.Int64Dtype(),
        'BAs' : float, 'NPs' : pd.Int64Dtype(),
        'BAa' : float, 'NPa' : pd.Int64Dtype(),
        'Rm' : float, 'RGm' : float, 'SRGm' : float,
        'Rs' : float, 'RGs' : float, 'SRGs' : float,
        'Ra' : float, 'RGa' : float, 'SRGa' : float,
        'CCSm' : float, 'CCPm' : float, 'ZCCPm' : float,
        'CCSs' : float, 'CCPs' : float, 'ZCCPs' : float,
        'CCSa' : float, 'CCPa' : float, 'ZCCPa' : float,
        'ZOm' : float, 'ZDm' : float, 'ZD-m' : float, 'ZD+m' : float,
        'ZOs' : float, 'ZDs' : float, 'ZD-s' : float, 'ZD+s' : float,
        'ZOa' : float, 'ZDa' : float, 'ZD-a' : float, 'ZD+a' : float,
    }

    def __init__(self,
        parse_logfile = False,
        ):

        self.parse_logfile = bool(parse_logfile)

    def __call__(self,
        dataset,
        f_label = None,
        ):

        scores_df, log_text = self.run_edstats_basic(
            dataset = dataset,
            f_label = f_label,
        )

        global_scores = (
            self._parse_log(
                log_text = log_text,
                )
            if (self.parse_logfile is True)
            else None
            )

        return EdstatsResults(
            results_table = scores_df,
            global_scores = global_scores,
        )

    def run_edstats_basic(self,
        dataset,
        f_label = None,
        ):

        import tempfile

        if not os.path.exists(dataset.model.filename):

            raise IOError(
                'PDB file for edstats does not exist! {!s}'.format(dataset.model.filename)
                )

        if not os.path.exists(dataset.data.filename):

            raise IOError(
                'MTZ file for edstats does not exist! {!s}'.format(dataset.data.filename)
                )

        # Get amplitude column
        if f_label is not None:
            self._check_valid_f_label(f_label, dataset)
        else:
            f_label = self._get_f_label(dataset)

        table_str_io, log_text = self._run_edstats_in_temp_dir_and_return_buffer(
            dataset = dataset,
            f_label = f_label,
            )

        df = self._extract_and_validate_table(table_str_io)

        return df, log_text

    def _check_valid_f_label(self, f_label, dataset):

        # Use column labels if given
        if (f_label is not None) and (f_label not in dataset.data.crystal.column_labels):
            raise Sorry(
                'Selected f_label ({}) not found in mtz file ({}) -- mtz contains columns {}'.format(
                    f_label,
                    dataset.data.filename,
                    dataset.data.crystal.column_labels,
                )
            )

    def _get_f_label(self, dataset):

        ###
        # change when MtzSummary code is updated (TODO)
        #
        from giant.io.mtz import MtzSummary
        #
        mtz_summary = MtzSummary(dataset.data.filename)
        #
        f_label = mtz_summary.labels.f
        #
        ###

        # Check for f_label
        if (not f_label):
            raise Sorry(
                'No F label selected/found in mtz file: {!s} -- mtz contains columns {}'.format(
                    dataset.data.filename,
                    dataset.data.crystal.column_labels,
                )
            )

        return f_label

    def _run_edstats_in_temp_dir_and_return_buffer(self, dataset, f_label):

        # Create a file handle and path for the output
        temp_handle, temp_path = tempfile.mkstemp(
            suffix = '.table',
            prefix = 'edstats_',
            )

        # Initialise Command Manager to run edstats
        from giant.dispatcher import (
            Dispatcher,
            )

        prog = Dispatcher('edstats.pl')

        prog.extend_args([
            '-hklin', dataset.data.filename,
            '-xyzin', dataset.model.filename,
            '-output', temp_path,
            '-flabel', f_label,
            '-noerror',
        ])

        # Log command
        logger(prog.as_string())

        # Run!
        prog.run()

        # Read the output
        with open(temp_path, 'r') as f:
            table_str = f.read()

        # Extract log_text for potential parsing
        log_text = str(prog.result.stdout)

        # Convert to buffer for reading
        import io
        table_str_io = io.StringIO(
            str(table_str)
            )

        # Cleanup
        os.remove(temp_path)

        return table_str_io, log_text

    def _extract_and_validate_table(self, buffer):

        # 0   RT:   Residue type (3-letter code).
        # 1   CI:   2 characters: chain ID (reset to underscore if blank) and alternate location indicator (if -usealt is specified).
        # 2   RN:   Residue name (including insertion code if present).
        # 39. MN:    Model number.
        # 40. CP:    PDB chain ID.
        # 41. NR:    PDB residue no.

        # Convert to dataframe
        df = pd.read_csv(
            buffer,
            #index_col = (39, 40, 2),  # MODEL, CHAIN, RESID
            delim_whitespace = True,
            dtype = self.EDSTATS_DTYPES,
        )

        # Assign index (must be done separately to prevent conversion to float)
        df = df.set_index(['MN', 'CP', 'RN'])

        assert list(df.index.names) == ['MN', 'CP', 'RN']

        assert list(df.columns) == (
            self.EDSTATS_COLUMNS[0:2] +
            self.EDSTATS_COLUMNS[3:39] +
            self.EDSTATS_COLUMNS[41:]
            )

        return df
