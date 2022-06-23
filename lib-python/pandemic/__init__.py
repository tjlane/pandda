import os, sys
import pandemic.resources
import giant

__version__ = "0.1.0"

HEADER_TEXT = """
------------------------------------------------------------------>
-
-  PanDEMIC Version {!s}""".format(__version__)+"""
-  Package Version {!s}""".format(giant.__version__)+"""
-
-  __________                    .___             .__
-  \\______   \\_____    ____    __| _/____   _____ |__| ____
-   |     ___/\\__  \\  /    \\  / __ |/ __ \\ /     \\|  |/ ___\\
-   |    |     / __ \\|   |  \\/ /_/ \\  ___/|  Y Y  \\  \\  \\___
-   |____|    (____  /___|  /\\____ |\\___  >__|_|  /__|\\___  >
-                  \\/     \\/      \\/    \\/      \\/        \\/
-
-  PAN-Dataset Ensemble Modelling of Iso-structural Crystals
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


class module_info(object):
    name        = 'pandemic'
    version     = __version__
    header_text = HEADER_TEXT

