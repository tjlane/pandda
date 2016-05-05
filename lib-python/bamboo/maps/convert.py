import os,sys

from bamboo.maps.utils import select_factors_and_phases_for_map
from bamboo.utils.mtz import get_mtz_summary_dict
from bamboo.wrappers.maputils import fft_mtz_to_map

def convert_mtz_to_map(mtzfile, mapfile=None, method='fft', maptype='2FOFC', force=False):
    """Converts an MTZ Format File to a MAP File (using fft as default)"""

    # Create mapfile if not given
    if not mapfile:
        mapfile = os.path.splitext(mtzfile)[0]+'.'+maptype+'.ccp4'
    # Check for file validity
    if os.path.exists(mapfile):
        if force:
            os.remove(mapfile)
        else:
            return None

    # Read the Column Headings from the MTZ File
    cols = get_mtz_summary_dict(mtzfile)['colheadings']
    selcols = select_factors_and_phases_for_map(cols,maptype)

    if method == 'fft':
        writer = fft_mtz_to_map(mtzfile, mapfile, selcols)
    else:
        raise Exception('Program not recognised: {!s}'.format(method))

    return writer
