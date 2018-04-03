import os, sys
import pandemic.resources

VERSION = '0.1.1'

HEADER_TEXT = """
------------------------------------------------------------------>
-
-  Pandemic Version {!s}""".format(VERSION)+"""
-
-  __________                    .___             .__
-  \\______   \\_____    ____    __| _/____   _____ |__| ____
-   |     ___/\\__  \\  /    \\  / __ |/ __ \\ /     \\|  |/ ___\\
-   |    |     / __ \\|   |  \\/ /_/ \\  ___/|  Y Y  \\  \\  \\___
-   |____|    (____  /___|  /\\____ |\\___  >__|_|  /__|\\___  >
-                  \\/     \\/      \\/    \\/      \\/        \\/
-
-  PAN-Dataset Ensemble Modelling of Isomorphous Crystals
-
------------------------------------------------------------------>

  When using this software, please cite:

       <insert exciting paper title here!>
       <insert exciting journal name here>

       N. Pearce & P. Gros (Unpublished).

  and include "PANDEMIC" in the KEYWRD sections of PDB depositions

------------------------------------------------------------------>

> {program}
{description}
------------------------------------------------------------------>
"""

LOGO_PATH = os.path.join(os.path.realpath(pandemic.resources.__path__[0]), 'pandemic-logo-small.png')

class module_info:
    name        = 'pandemic'
    version     = VERSION
    header_text = HEADER_TEXT

