import os, sys
import pkg_resources

__version__ = '0.2.12'

from libtbx import group_args, adopt_init_args

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

