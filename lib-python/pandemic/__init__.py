import os, sys
import pandemic.resources

import bamboo
__version__ = bamboo.__version__

HEADER_TEXT = """
------------------------------------------------------------------>
-
-  Package Version {!s}""".format(__version__)+"""
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
    version     = __version__
    header_text = HEADER_TEXT

