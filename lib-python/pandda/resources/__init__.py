import os

_this_path = os.path.realpath(__path__[0])

FULL_LOGO_PATH = os.path.join(_this_path, 'pandda-logo.png')
SMALL_LOGO_PATH = os.path.join(_this_path, 'pandda-logo-small.png')

# LOGO_TEXT = """
# __________                    .___  .___
# \\______   \\_____    ____    __| _/__| _/____   
#  |     ___/\\__  \\  /    \\  / __ |/ __ |\\__  \\ 
#  |    |     / __ \\|   |  \\/ /_/ / /_/ | / __ \\_
#  |____|    (____  /___|  /\\____ \\____ |(____  /
#                 \\/     \\/      \\/    \\/     \\/ 
# """

FORMAT_LOGO_TEXT = """
__________                    .___  .___
\\______   \\_____    ____    __| _/__| _/____   
 |     ___/\\__  \\  /    \\  / __ |/ __ |\\__  \\ 
 |    |     / __ \\|   |  \\/ /_/ / /_/ | / __ \\_
 |____|    (____  /___|  /\\____ \\____ |(____  / {}
                \\/     \\/      \\/    \\/     \\/ 
"""

LOGO_TEXT = FORMAT_LOGO_TEXT.format('')