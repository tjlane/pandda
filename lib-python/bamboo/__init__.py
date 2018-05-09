# Contains all non-specialist code for the analysis of protein and ligand structures
# Does NOT require cctbx

import os, sys

from bamboo.common.errors import *

try:
    VERSION = pkg_resources.get_distribution("panddas").version
except:
    VERSION = '(developer -- see setup.py file)'
