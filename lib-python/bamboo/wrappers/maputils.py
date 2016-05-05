import os,sys

from bamboo.common.command import CommandManager

def fft_mtz_to_map(mtz_file, map_file, cols):
    """Converts an MTZ Format File to a MAP File (using fft as default)"""

    # Initialise
    writer = CommandManager('fft')
    # Set Program Arguments
    writer.add_command_line_arguments('hklin',mtz_file,'mapout',map_file)
    # Set Program Input
    writer.add_standard_input(['LABIN F1={!s} PHI={!s}'.format(cols['F'],cols['P']),'END'])
    # RUN!
    writer.run()
    # Check Output
    if writer.process.returncode != 0:
        print('\nOUT\n\n'+writer.out)
        print('\nERR\n\n'+writer.err)
        raise RuntimeError('fft failed to generate map from {!s}'.format(mtz_file))

    return writer

def mask_map(mapin, maskpdb, mapout, border=1):
    """Takes mapin and masks around atoms in maskpdb"""

    # Masking object
    masker = CommandManager('mapmask')
    # Set input files
    masker.add_command_line_arguments(['mapin',mapin,'mapout',mapout,'xyzin',maskpdb])
    # Set stdin
    masker.add_standard_input(['BORDER {!s}'.format(border),'END'])
    # Run!
    masker.run()
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
    masker.add_command_line_arguments(['mapin',mapin,'mapout',mapout])
    # Set stdin
    masker.add_standard_input(['XYZLIM ASU','END'])
    # Run!
    masker.run()
    # Report errors
    if masker.process.returncode!=0:
        raise RuntimeError('mapmask failed to create asu map from {!s}'.format(mapin))

    # Return Command Managers for flexible handling of out & err
    return masker

