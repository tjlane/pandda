import os, sys
import pkg_resources

PANDDA_LIB_TOP = os.path.realpath(os.path.abspath(os.path.join(__path__[0], '..')))
try:
    PANDDA_VERSION = pkg_resources.get_distribution("panddas").version
except:
    PANDDA_VERSION = 'xxx'

PANDDA_TEXT = """
------------------------------------------------------------------>
-
-  Pandda Version {!s}
-
-  __________                    .___  .___
-  \\______   \\_____    ____    __| _/__| _/____    ______
-   |     ___/\\__  \\  /    \\  / __ |/ __ |\\__  \\  /  ___/
-   |    |     / __ \\|   |  \\/ /_/ / /_/ | / __ \\_\\___ \\
-   |____|    (____  /___|  /\\____ \\____ |(____  /____  >
-                  \\/     \\/      \\/    \\/     \\/     \\/
-
-  PAN-Dataset Density Analysis of Crystallographic Fragment Screens
-
-  When using this software, please cite:
-
-       Unpublished,
-       N Pearce, et al.
-
-  and include "PANDDA" in the KEYWRD section of related PDB files
-
------------------------------------------------------------------>
""".format(PANDDA_VERSION)

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
        print 'Hi {!s}. Welcome to Pandda.'.format(user.upper())
    except: pass

