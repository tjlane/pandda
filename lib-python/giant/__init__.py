import os, sys
import pkg_resources

import bamboo
__version__ = bamboo.__version__

HEADER_TEXT = """
------------------------------------------------------------------>
-
-  Package Version {!s}""".format(__version__)+"""
-
-         .__               __
-    ____ |__|____    _____/  |_
-   / ___\\|  \\__  \\  /    \\   __\\
-  / /_/  >  |/ __ \\|   |  \\  |
-  \\___  /|__(____  /___|  /__|     A crystallographic toolbox for
-  /_____/         \\/     \\/        protein structure determination
-
------------------------------------------------------------------>

> {program}
{description}
------------------------------------------------------------------>
"""

class module_info:
    name        = 'giant'
    version     = __version__
    header_text = HEADER_TEXT
