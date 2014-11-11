#! /usr/local/python/python2.7.3-64bit/bin/python

import os, sys, shutil

from Bamboo.Utils.MTZ import get_mtz_summary
from Bamboo.Wrappers.PdbUtils import create_cryst_line

def create_cryst_line_from_mtz(pdbin, pdbout, mtz):
    """Takes information from the mtz file and creates a cryst line in pdbin"""

    # Get the spacegroup number from the mtzfile
    mtzsummary = get_mtz_summary(mtz)
    sg   = mtzsummary['spacegroup'].replace(' ','')
    cell = mtzsummary['cell']

    pdbset = create_cryst_line(pdbin, pdbout, sg, cell)

    return pdbset
