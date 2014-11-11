#! /usr/local/python/python2.7.3-64bit/bin/python

# NEW Python Library for the PANDDA project

# Contains all non-specialist code for the analysis of protein and ligand structures

import os, sys

if '/usr/local/scripts/' not in sys.path:
    sys.path.append('/usr/local/scripts')
if '/home/npearce/bin/PANDDAs' not in sys.path:
    sys.path.append('/home/npearce/bin/PANDDAs')
if '/home/npearce/bin/' not in sys.path:
    sys.path.append('/home/npearce/bin')

from Bamboo.Common.Errors import *


