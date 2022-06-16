try:
    import giant.logs as lg
except ImportError:
    import logging as lg
logger = lg.getLogger(__name__)

import os, sys
import pkg_resources

#from libtbx import group_args, adopt_init_args

__version_tuple__ = (1, 0, 0)
__version__ = '.'.join(map(str,__version_tuple__))

HEADER_TEXT = """
------------------------------------------------------------------>
-
-  Package Version {version}
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

class module_info(object): # backwards compatibility
    name = 'giant'
    version = __version__
    header_text = HEADER_TEXT.format(
        program = "{program}",
        description = "{description}",
        version = __version__,
        )


class ModuleBanner(object):

    banner_string = HEADER_TEXT.format(
        version = "{module_version}",
        program = "{program}",
        description = "{description}",
        )

    module_name = "giant"
    module_version = __version__

    def __init__(self,
        program = None,
        description = None,
        ):

        if (program is None):
            program = self.module_name

        if (description is None):
            description = ""

        self.program = program
        self.description = description

    def __str__(self):

        s = self.banner_string.format(
            module_name = self.module_name,
            module_version = self.module_version,
            program = self.program,
            description = self.description,
            )

        return s

module_banner = ModuleBanner()
