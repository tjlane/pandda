import giant.logs as lg
logger = lg.getLogger(__name__)

import os, sys
import pandas

from giant.exceptions import Sorry, Failure

from libtbx import adopt_init_args


class EdstatsFactory:

    # Output columns in order
    edstats_columns = [
        'RT', 'CI', 'RN',
        'BAm', 'NPm', 'Rm', 'RGm', 'SRGm', 'CCSm', 'CCPm', 'ZCCPm', 'ZOm', 'ZDm', 'ZD-m', 'ZD+m',
        'BAs', 'NPs', 'Rs', 'RGs', 'SRGs', 'CCSs', 'CCPs', 'ZCCPs', 'ZOs', 'ZDs', 'ZD-s', 'ZD+s',
        'BAa', 'NPa', 'Ra', 'RGa', 'SRGa', 'CCSa', 'CCPa', 'ZCCPa', 'ZOa', 'ZDa', 'ZD-a', 'ZD+a',
        'MN', 'CP', 'NR'
    ]

    # Output datatypes
    edstats_dtypes = {
        'RT' : str, # Keep as string for use as labels
        'CI' : str, # Keep as string for use as labels
        'RN' : str, # Keep as string for use as labels
        'MN' : str, # Keep as string for use as labels
        'CP' : str, # Keep as string for use as labels
        'NR' : str, # Keep as string for use as labels
        'BAm' : float, 'NPm' : pandas.Int64Dtype(),
        'BAs' : float, 'NPs' : pandas.Int64Dtype(),
        'BAa' : float, 'NPa' : pandas.Int64Dtype(),
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
        adopt_init_args(self, locals())

    def __call__(self,
        pdb_file,
        mtz_file,
        f_label = None,
        ):

        scores_df, dispatcher = self.run_edstats_basic(
            pdb_file = pdb_file,
            mtz_file = mtz_file,
            f_label = f_label,
        )

        log_data = None
        if (self.parse_logfile is True):
            log_data = EdstatsLogProcessor(logtext=dispatcher.result.stdout)

        return EdstatsResults(
            scores_df = scores_df,
            log_data = log_data,
        )

    def run_edstats_basic(self,
        pdb_file,
        mtz_file,
        f_label = None,
        ):

        import tempfile

        if not os.path.exists(pdb_file):
            raise IOError('PDB file for edstats does not exist! {!s}'.format(pdb_file))
        if not os.path.exists(mtz_file):
            raise IOError('MTZ file for edstats does not exist! {!s}'.format(mtz_file))

        # Create a file handle and path for the output
        temp_handle, temp_path = tempfile.mkstemp(suffix='.table', prefix='edstats_')

        # Collate summary of MTZ file
        from giant.io.mtz import MtzSummary
        m_summ = MtzSummary(mtz_file)

        # Use column labels if given
        if (f_label is not None) and (f_label not in m_summ.summary_dict['colheadings']):
            raise Sorry(
                'Selected f_label ({}) not found in mtz file ({}) -- mtz contains columns {}'.format(
                    f_label,
                    mtz_file,
                    m_summ.summary_dict['colheadings'],
                )
            )
        # else guess the labels in the mtzfile
        else:
            f_label = m_summ.labels.f

        # Check for f_label
        if (not f_label):
            raise Sorry(
                'No F label selected/found in mtz file: {!s} -- mtz contains columns {}'.format(
                    mtz_file,
                    m_summ.summary_dict['colheadings'],
                )
            )

        # Initialise Command Manager to run edstats
        from giant.dispatcher import Dispatcher
        prog = Dispatcher('edstats.pl')
        prog.extend_args([
            '-hklin', mtz_file,
            '-xyzin', pdb_file,
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
        os.remove(temp_path)

        # Convert to buffer for reading
        import io
        table_str_io = io.StringIO(unicode(table_str))

        # 0   RT:   Residue type (3-letter code).
        # 1   CI:   2 characters: chain ID (reset to underscore if blank) and alternate location indicator (if -usealt is specified).
        # 2   RN:   Residue name (including insertion code if present).
        # 39. MN:    Model number.
        # 40. CP:    PDB chain ID.
        # 41. NR:    PDB residue no.

        # Convert to dataframe
        df = pandas.read_csv(
            table_str_io,
            #index_col = (39, 40, 2),  # MODEL, CHAIN, RESID
            delim_whitespace = True,
            dtype = self.edstats_dtypes,
        )

        # Assign index (must be done separately to prevent conversion to float)
        df = df.set_index(['MN', 'CP', 'RN'])

        assert list(df.index.names) == ['MN', 'CP', 'RN']
        assert list(df.columns) == (self.edstats_columns[0:2] + self.edstats_columns[3:39] + self.edstats_columns[41:])

        return df, prog


class EdstatsResults:

    def __init__(self,
        scores_df,
        log_data = None,
        ):

        self.scores_df = scores_df
        self.log_data = log_data

    def get_label(self, residue_group):
        rg = residue_group
        model_no = rg.parent().parent().id.strip()
        chain_id = rg.parent().id.strip()
        res_id = rg.resid().strip()
        if model_no == '':
            model_no = '1'
        if chain_id == '':
            chain_id = '_'
        return (model_no, chain_id, res_id)

    def extract_residue_group_scores(self, residue_group=None, residue_group_label=None):
        """Extract density quality metrics for a residue group from precalculated edstats scores"""

        if [residue_group, residue_group_label].count(None) != 1:
            raise ValueError('Provide EITHER residue_group OR residue_group_label')

        # Shortcuts
        rg = residue_group
        rg_label = residue_group_label

        if (rg_label is None):
            rg_label = self.get_label(residue_group=rg)

        # Check label associated with a record
        if not self.scores_df.index.contains(rg_label):
            raise KeyError('The given label does not have a corresponding result: {}'.format(rg_label))

        # Extract residue scores
        ed_scores = self.scores_df.loc[rg_label]

        # Append scores to data_table
        output_dict = {
            'rscc' : ed_scores['CCSa'],
            'rsr'  : ed_scores['Ra'],
            'b_av' : ed_scores['BAa'],
            'rszo' : ed_scores['ZOa'],
            'rszd' : ed_scores['ZDa'],
        }

        return output_dict


class EdstatsLogProcessor:
    """Class to process, store and present the output from edstats.pl"""

    def __init__(self, logtext):

        self.logtext = logtext

        self.mean_biso = None

        self.rms_scores = {
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

        # Populate fields
        self._parselogtext()

    def _parselogtext(self):
        """Read and sort through log input"""

        output = [line for line in self.logtext.split('\n') if line]
        for i, line in enumerate(output):
            if line.startswith('Overall mean Biso:'):
                self.mean_biso = float(line.replace('Overall mean Biso:',''))
            elif line.startswith('RMS Z-scores for protein residues:'):
                assert ['Main','Side','All'] == output[i+1].strip().split(), 'COLUMN HEADINGS ARE NOT WHAT IS EXPECTED! {!s}'.format(line)
                for incr in [2,3,4]:
                    score, main_val, side_val, all_val = output[i+incr].strip().split()
                    self.rms_scores['_'.join(['atm','main',score])] = float(main_val)
                    self.rms_scores['_'.join(['atm','side',score])] = float(side_val)
                    self.rms_scores['_'.join(['atm','all',score])] = float(all_val)
            elif line.startswith('RMS Z-scores for hetero residues:'):
                assert ['Main','Side','All'] == output[i+1].strip().split(), 'COLUMN HEADINGS ARE NOT WHAT IS EXPECTED! {!s}'.format(line)
                for incr in [2,3,4]:
                    score, main_val, side_val, all_val = output[i+incr].strip().split()
                    self.rms_scores['_'.join(['het','main',score])] = float(main_val)
                    self.rms_scores['_'.join(['het','side',score])] = float(side_val)
                    self.rms_scores['_'.join(['het','all',score])] = float(all_val)

########################################################################################################

