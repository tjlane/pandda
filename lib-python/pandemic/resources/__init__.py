import os

_this_path = os.path.realpath(__path__[0])

FULL_LOGO_PATH = os.path.join(_this_path, 'pandemic-logo.png')
SMALL_LOGO_PATH = os.path.join(_this_path, 'pandemic-logo-small.png')
NO_IMAGE_PATH = os.path.join(_this_path, 'pandemic-no-image.png')
NO_IMAGE_PATH_ADP = os.path.join(_this_path, 'pandemic-adp-no-image.png')

# LOGO_TEXT = """
# __________                     ___              __
# \\______   \\_____    ____    __| _/____   _____ |__| ____
#  |     ___/\\__  \\  /    \\  / __ |/ __ \\ /     \\|  |/ ___\\
#  |    |     / __ \\|   |  \\/ /_/ \\  ___/|  Y Y  \\  \\  \\___
#  |____|    (____  /___|  /\\____ |\\___  >__|_|  /__|\\___  >
#                 \\/     \\/      \\/    \\/      \\/        \\/
# """

FORMAT_LOGO_TEXT = """
__________                     ___              __
\\______   \\_____    ____    __| _/____   _____ |__| ____
 |     ___/\\__  \\  /    \\  / __ |/ __ \\ /     \\|  |/ ___\\
 |    |     / __ \\|   |  \\/ /_/ \\  ___/|  Y Y  \\  \\  \\___
 |____|    (____  /___|  /\\____ |\\___  >__|_|  /__|\\___  > {}
                \\/     \\/      \\/    \\/      \\/        \\/
"""

LOGO_TEXT = FORMAT_LOGO_TEXT.format('')