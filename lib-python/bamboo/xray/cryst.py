import os, sys, shutil

from bamboo.utils.mtz import get_mtz_summary_dict
from bamboo.wrappers.pdbutils import create_cryst_line

def create_cryst_line_from_mtz(pdbin, pdbout, mtz):
    """Takes information from the mtz file and creates a cryst line in pdbin"""

    # Get the spacegroup number from the mtzfile
    mtzsummary = get_mtz_summary_dict(mtz)
    sg   = mtzsummary['spacegroup'].replace(' ','')
    cell = mtzsummary['cell']

    pdbset = create_cryst_line(pdbin, pdbout, sg, cell)

    return pdbset
