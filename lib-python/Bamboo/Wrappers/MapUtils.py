#! /usr/local/python/python2.7.3-64bit/bin/python

import os,sys

from Bamboo.Common.Command import CommandManager

def fft_mtz_to_map(mtzfile, mapfile, cols):
    """Converts an MTZ Format File to a MAP File (using fft as default)"""

    # Initialise
    writer = CommandManager('fft')
    # Set Program Arguments
    writer.SetArguments('hklin',mtzfile,'mapout',mapfile)
    # Set Program Input
    writer.SetInput(['LABIN F1={!s} PHI={!s}'.format(cols['F'],cols['P']),'END'])
    # RUN!
    writer.Run()
    # Check Output
    if writer.process.returncode != 0:
        print('\nOUT\n\n'+writer.out)
        print('\nERR\n\n'+writer.err)
        raise RuntimeError('fft failed to generate map from {!s}'.format(mtzfile))

    return writer

def mask_map(mapin, maskpdb, mapout, border=1):
    """Takes mapin and masks around atoms in maskpdb"""

    # Masking object
    masker = CommandManager('mapmask')
    # Set input files
    masker.SetArguments(['mapin',mapin,'mapout',mapout,'xyzin',maskpdb])
    # Set stdin
    masker.SetInput(['BORDER {!s}'.format(border),'END'])
    # Run!
    masker.Run()
    # Report errors
    if masker.process.returncode!=0:
        raise RuntimeError('mapmask failed to mask map {!s}'.format(mapin))

    # Return Command Managers for flexible handling of out & err
    return masker

def create_asu_map(mapin, mapout):
    """Takes mapin and masks to the asymmetric unit"""

    # Masking object
    masker = CommandManager('mapmask')
    # Set input files
    masker.SetArguments(['mapin',mapin,'mapout',mapout])
    # Set stdin
    masker.SetInput(['XYZLIM ASU','END'])
    # Run!
    masker.Run()
    # Report errors
    if masker.process.returncode!=0:
        raise RuntimeError('mapmask failed to create asu map from {!s}'.format(mapin))

    # Return Command Managers for flexible handling of out & err
    return masker

