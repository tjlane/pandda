import giant.logs as lg
logger = lg.getLogger(__name__)

import os, re, shutil, tempfile

from libtbx import adopt_init_args

from giant.dispatcher import Dispatcher
from giant.exceptions import Failure

# 2FOFC cols
F_COMP_OPTIONS = ['2FOFCWT','FWT']
PHI_COMP_OPTIONS = ['PH2FOFCWT','PHWT','PHFWT']

# FOFC cols
F_DIFF_OPTIONS = ['FOFCWT','DELFWT']
PHI_DIFF_OPTIONS = ['PHFOFCWT','DELPHWT','PHDELWT']

def _get_first_if_not_none(item):
    if (item is not None):
        if len(item) == 0:
            return None
        return item[0]
    return item

def _get_regex_matches(text, regex_str, n=1):
    import re
    regex = re.compile(regex_str)
    matches = regex.findall(text)
    if (n is not None):
        assert len(matches) == n
    if n == 1:
        return matches[0]
    return matches

def assign_column_labels(col_labs, col_types):
    """Attempt to automatically assign meaning to input column names"""

    lab_dict = {}

    # Sort the flags
    assert len(col_labs) == len(col_types)
    # F: Amplitudes, H: Miller Indices, I: Integers, P: Phases, Q: Standard Deviations, J: Intensities
    sfac = [col_labs[i] for i,t in enumerate(col_types) if t=='F']
    inty = [col_labs[i] for i,t in enumerate(col_types) if t=='J']
    mill = [col_labs[i] for i,t in enumerate(col_types) if t=='H']
    ints = [col_labs[i] for i,t in enumerate(col_types) if t=='I']
    phas = [col_labs[i] for i,t in enumerate(col_types) if t=='P']
    sdev = [col_labs[i] for i,t in enumerate(col_types) if t=='Q']
    foms = [col_labs[i] for i,t in enumerate(col_types) if t=='W']

    # Record 2FOFC pair (composite map)
    lab_dict['f_comp_labels'] = [s for s in sfac if s in F_COMP_OPTIONS]
    lab_dict['phi_comp_labels'] = [p for p in phas if p in PHI_COMP_OPTIONS]
    # Record FOFC pair (difference map)
    lab_dict['f_diff_labels'] = [s for s in sfac if s in F_DIFF_OPTIONS]
    lab_dict['phi_diff_labels'] = [p for p in phas if p in PHI_DIFF_OPTIONS]

    # Find the main F col
    lab_dict['f_labels'] = [
        s for s in sfac if (
            s in ['F','FP','FCTR','FOSC','FOBS']
            or s.startswith('F_')
            or s.startswith('FP_')
            or s.upper().startswith('F-OBS')
            or s.upper().startswith('FOUT_')
        )
    ]
    # Find the main SIGF col
    lab_dict['sigf_labels'] = [
        s for s in sdev if (
            s in ['SIGF','SIGFP','SIGFOBS']
            or s.startswith('SIGF_')
            or s.startswith('SIGFP_')
            or s.upper().startswith('SIGF-OBS')
            or s.upper().startswith('SIGFOUT_')
        )
    ]
    # Find the I cols
    lab_dict['i_labels'] = [
        s for s in inty if (
            s.startswith('I')
        )
    ]
    # Find the SIGI cols
    lab_dict['sigi_labels'] = [
        s for s in sdev if (
            s.startswith('SIGI')
        )
    ]
    # Find the F_calc
    lab_dict['f_calc_labels'] = [
        s for s in sfac if (
            s in ['FC', 'FMODEL']
            or s.startswith('FC_')
            or s.upper().startswith('F-MODEL')
        )
    ]
    # Find the PHI_calcs
    lab_dict['phi_calc_labels'] = [
        p for p in phas if (
            p in ['PHIC', 'PHIFC', 'PHIFMODEL']
            or p.startswith('PHIC_')
            or p.upper().startswith('PHIF-MODEL')
        )
    ]
    # Find the RFree Flag
    lab_dict['r_free_labels'] = [
        r for r in ints if (
            'FREE' in r.upper()
        )
    ]
    # Figure-o-merits
    lab_dict['fom_labels'] = [
        l for l in foms if (
            'FOM' in l.upper()
        )
    ]

    used_labels = []; [used_labels.extend(vals) for vals in list(lab_dict.values())]
    used_labels = sorted(used_labels, key=lambda x: col_labs.index(x))
    unused_labels = [l for l in col_labs if (l not in used_labels)]

    lab_dict['assigned_labels'] = used_labels
    lab_dict['unassigned_labels'] = unused_labels

    return lab_dict

class MtzHeaderData(object):
    """Object to hold the meta data for an MTZ file"""

    def __init__(self,
        resolution_low = None,
        resolution_high = None,
        spacegroup = None,
        spacegroup_no = None,
        cell = None,
        ):

        cell = tuple(cell)

        adopt_init_args(self, locals())

    def __str__(self):

        s = """
        Resolution Range: {low} - {high} A
        Spacegroup: {sg} (No. {sgno})
        Cell: {cell}
        """.format(
            low = self.resolution_low,
            high = self.resolution_high,
            sg = self.spacegroup,
            sgno = self.spacegroup_no,
            cell = self.cell,
        )
        return s

    def summary(self):
        return str(self)

    def as_cryst(self):
        return 'CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f}{:>11}{:4}'.format(
            self.cell[0],
            self.cell[1],
            self.cell[2],
            self.cell[3],
            self.cell[4],
            self.cell[5],
            self.spacegroup,
            ' ',
        )


class MtzColumnLabels(object):
    """Object to contain the column labels of an MTZ file"""

    def __init__(self,
        f = None,
        sigf = None,
        i = None,
        sigi = None,
        free = None,
        fom = None,
        f_calc = None,
        phi_calc = None,
        f_comp = None,
        phi_comp = None,
        f_diff = None,
        phi_diff = None,
        ):

        adopt_init_args(self, locals())

    def __str__(self):
        s = """
        MTZ Column labels:
            F: {f} ({sigf})
            I: {i} ({sigi})
            FREE: {free}
            FOM: {fom}
            CALC: {f_calc} / {phi_calc}
            COMP: {f_comp} / {phi_comp}
            DIFF: {f_diff} / {phi_diff}
        """.format(
            f = self.f,
            sigf = self.sigf,
            i = self.i,
            sigi = self.sigi,
            free = self.free,
            fom = self.fom,
            f_calc = self.f_calc,
            phi_calc = self.phi_calc,
            f_comp = self.f_comp,
            phi_comp = self.phi_comp,
            f_diff = self.f_diff,
            phi_diff = self.phi_diff,
        )
        return s

    def summary(self):
        return str(self)


class MtzSummary(object):
    """Class for summarising an MTZ file"""

    def __init__(self, mtz_file):

        summary_dict = get_mtz_summary_dict(mtz_file)

        self.header = MtzHeaderData(
            resolution_low  = summary_dict.get('reslow'),
            resolution_high = summary_dict.get('reshigh'),
            spacegroup      = summary_dict.get('spacegroup'),
            spacegroup_no   = summary_dict.get('spacegroupno'),
            cell            = summary_dict.get('cell'),
        )

        # Sort column_labels
        col_labs_dict = assign_column_labels(
            col_labs = summary_dict.get('colheadings'),
            col_types = summary_dict.get('coltypes'),
        )

        self.labels = MtzColumnLabels(
            f           = _get_first_if_not_none(col_labs_dict.get('f_labels')),
            sigf        = _get_first_if_not_none(col_labs_dict.get('sigf_labels')),
            i           = _get_first_if_not_none(col_labs_dict.get('i_labels')),
            sigi        = _get_first_if_not_none(col_labs_dict.get('sigi_labels')),
            free        = _get_first_if_not_none(col_labs_dict.get('r_free_labels')),
            fom         = _get_first_if_not_none(col_labs_dict.get('fom_labels')),
            f_calc      = _get_first_if_not_none(col_labs_dict.get('f_calc_labels')),
            phi_calc    = _get_first_if_not_none(col_labs_dict.get('phi_calc_labels')),
            f_comp      = _get_first_if_not_none(col_labs_dict.get('f_comp_labels')),
            phi_comp    = _get_first_if_not_none(col_labs_dict.get('phi_comp_labels')),
            f_diff      = _get_first_if_not_none(col_labs_dict.get('f_diff_labels')),
            phi_diff    = _get_first_if_not_none(col_labs_dict.get('phi_diff_labels')),
        )

        self.summary_dict = summary_dict

    def __str__(self):

        s = [
            self.header.summary(),
            self.labels.summary(),
        ]

        return '\n'.join(s)

    def summary(self):
        return str(self)


def get_mtz_summary_dict(mtz_file):
    """Get fields from an MTZ"""

    # Extract Contents of MTZ
    mtzdmp = Dispatcher('mtzdmp')
    mtzdmp.append_arg(mtz_file)
    mtzdmp.run()

    # Check for errors
    if mtzdmp.result['exitcode'] != 0:
        raise RuntimeError(
            'mtzdmp failed to read file {!s}\n> STDOUT \n{!s}\n> STDERR\n{!s}'.format(
                mtz_file,
                str(mtzdmp.result.stdout),
                str(mtzdmp.result.stderr),
            )
        )

    return parse_mtzdmp_logtext(logtext=str(mtzdmp.result.stdout))

def parse_mtzdmp_logtext(logtext):

    # Create empty dict to contain the summary
    summary = {}

    # Get the resolution range
    res_range_str = _get_regex_matches(
        text = logtext,
        regex_str = r'\*  Resolution Range :.*\\n.*\\n.*\((.*)A \)\\n',
    )
    res_range = list(map(float, res_range_str.replace(' ','').split('-')))
    summary['reslow'], summary['reshigh'] = res_range

    # Get the Number of Columns
    n_cols_str = _get_regex_matches(
        text = logtext,
        regex_str = r'\* Number of Columns =(.*?)\\n',
    )
    n_cols = int(n_cols_str.strip())
    summary['numcols'] = n_cols

    # Get the Number of Reflections
    n_refl_str = _get_regex_matches(
        text = logtext,
        regex_str = r'\* Number of Reflections =(.*?)\\n',
    )
    n_refl = int(n_refl_str.strip())
    summary['numreflections'] = n_refl

    # Get the Column Labels
    col_labs_str = _get_regex_matches(
        text = logtext,
        regex_str = r'\* Column Labels :.*?\\n.*?\\n(.*?)\\n',
    )
    col_labs = col_labs_str.strip().split()
    summary['colheadings'] = col_labs

    # Get the Column Types
    col_types_str = _get_regex_matches(
        text = logtext,
        regex_str = r'\* Column Types :.*?\\n.*?\\n(.*?)\\n',
    )
    col_types = col_types_str.strip().split()
    summary['coltypes'] = col_types

    # Get the different datasets
    col_datasets_str = _get_regex_matches(
        text = logtext,
        regex_str = r'\* Associated datasets :.*?\\n.*?\\n(.*?)\\n',
    )
    col_datasets = col_datasets_str.strip().split()
    summary['coldatasets'] = col_datasets

    # Get the Spacegroup
    spacegroup_strs = _get_regex_matches(
        text = logtext,
        regex_str = r'\* Space group =.*?\\\'(.*?)\\\'.\(number(.*?)\)',
    )
    spacegroup = spacegroup_strs[0].strip()
    spacegroup_no = int(spacegroup_strs[1].strip())
    summary['spacegroup'] = spacegroup
    summary['spacegroupno'] = spacegroup_no

    # Get the Cell Dimensions
    cell_str = _get_regex_matches(
        text = logtext,
        regex_str = r'\* Cell Dimensions :.*?\\n.*?\\n(.*?)\\n',
    )
    cell = list(map(float,cell_str.split()))
    summary['cell'] = cell

    return summary

