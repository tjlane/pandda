import os, sys
import pkg_resources
import pandda.resources

try:
    VERSION = pkg_resources.get_distribution("panddas").version
except:
    VERSION = 'developer -- see setup.py file'

HEADER_TEXT = """
------------------------------------------------------------------>
-
-  Pandda Version {!s}""".format(VERSION)+"""
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

> {program}
{description}
------------------------------------------------------------------>
"""

LOGO_PATH = os.path.join(os.path.realpath(pandda.resources.__path__[0]), 'pandda-logo-small.png')

class module_info:
    name        = 'pandda'
    version     = VERSION
    header_text = HEADER_TEXT

def welcome(user=None):
    """Welcome message"""

    if user is None:
        try:    user = os.getlogin()
        except: user = "Unknown user"

    try:
        from greetings import get_greeting
        print get_greeting(username=user)
    except SystemExit: raise
    except: pass

    try:
        print '\nHi {!s}. Welcome to Pandda.'.format(user.upper())
    except: pass

