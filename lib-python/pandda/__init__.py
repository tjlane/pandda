import os, sys
import giant
import pandda.resources

__version_tuple__ = (2, 0, 0)
__version__ = '.'.join(map(str,__version_tuple__))

HEADER_TEXT = """
------------------------------------------------------------------>
-
-  Giant Package Version: {giant_version}
-  Pandda Version: {module_version}
-
-  __________                    .___  .___
-  \\______   \\_____    ____    __| _/__| _/____    ______
-   |     ___/\\__  \\  /    \\  / __ |/ __ |\\__  \\  /  ___/
-   |    |     / __ \\|   |  \\/ /_/ / /_/ | / __ \\_\\___ \\
-   |____|    (____  /___|  /\\____ \\____ |(____  /____  >
-                  \\/     \\/      \\/    \\/     \\/     \\/
-
-  PAN-Dataset Density Analysis of Crystallographic Datasets
-
------------------------------------------------------------------>

  When using this software, please cite:

       A Multi-Crystal Method for Extracting Obscured Signal
       from Crystallographic Electron Density.

       N Pearce, et al. Nature Communications (2017).

  and include "PANDDA" in the KEYWRD sections of PDB depositions

------------------------------------------------------------------>

> Program: {program}
{description}
------------------------------------------------------------------>
""".format(
    giant_version = str(giant.__version__),
    module_version = "{module_version}",
    program = "{program}",
    description = "{description}",
    )


class ModuleBanner(giant.ModuleBanner):
    
    banner_string = HEADER_TEXT

    module_name = 'pandda'
    module_version = __version__


module_banner = ModuleBanner()
