import os

from bamboo.common.command import CommandManager
from bamboo.maps.convert import convert_mtz_to_map
from bamboo.wrappers.maputils import mask_map

def create_masked_map_from_file_objects(mtzobj,pdbobj,mapout=None,maptype='2FOFC',force=False):
    """Takes the ligand described by ligobj and masks against the density from pdbobj & mtzobj"""

    # oriented and positioned structure of ligand only
    pdb_file = pdbobj.path
    # reflection data path
    mtz_file = mtzobj.path
    # Create the map
    fft, masker = create_masked_map(mtzin=mtz_file, maskpdb=pdb_file, mapout=mapout, maptype=maptype)

    return fft, masker

def create_masked_map(mtzin, maskpdb, mapout=None, maptype='2FOFC', border=2):
    """Creates density from mtzin, and masks around maskpdb - NOTE THAT MASKPDB MUST BE ALIGNED WITH MTZIN"""

    # Create mapfile if not given
    if not mapout:
        mapout = os.path.splitext(maskpdb)[0]+'.masked.'+maptype+'.ccp4'
    # Temporary (unmasked) map
    tempmap = mapout+'.nomask'
    # Convert the inital MAT to a MAP file
    fft = convert_mtz_to_map(mtzfile=mtzin, mapfile=tempmap, maptype=maptype)
    # Mask the map
    masker = mask_map(mapin=tempmap,maskpdb=maskpdb,mapout=mapout,border=border)
    # Remove the temporary mapfile
    os.remove(tempmap)
    # Return Command Managers for flexible handling of out & err
    return fft, masker


