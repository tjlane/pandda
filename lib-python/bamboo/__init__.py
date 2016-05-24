# Contains all non-specialist code for the analysis of protein and ligand structures
# Does NOT require cctbx

import os, sys

from bamboo.common.errors import *

def div_bar(width=30):
    if width < 5: width=5
    return '#'*2 + '='*(width-5) + '>'*3

