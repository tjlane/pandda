import os, sys

PANDEMIC_TOP = os.path.join(__path__[0], '..', '..')
PANDEMIC_VERSION = 0.0.0
PANDEMIC_TEXT = """
------------------------------------------------------------------>
-
-  Pandemic Version {!s}
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
-  When using this software, please cite:
-
-       Unpublished,
-       N Pearce, et al.
-
-  and include "PANDEMIC" in the KEYWRD section of related PDB files
-
------------------------------------------------------------------>
""".format(PANDEMIC_VERSION)

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
        print 'Hi {!s}. Welcome to Pandemic.'.format(user.upper())
    except: pass
