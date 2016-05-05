import os, sys

from bamboo.common.command import CommandManager
from bamboo.utils.mtz import MtzSummary

def convert_intensities_to_amplitudes(mtzin, mtzout):
    """Takes an input mtz and converts the intensities to structure factors for model building"""

    mtzobj = MtzSummary(mtzin)

    # Copy data columns from new mtz file
    I, SIGI = mtzobj.label.i, mtzobj.label.sigi
    if not (I and SIGI):
        raise LabelError('No Intensities found in {!s}'.format(mtzin))

    # TODO If rfree is present, retain rfree from the reference mtz
#    RFree = mtzobj.label.free

    # Initialise Commander
    CTRUNC = CommandManager('ctruncate')
    # Set command arguments
    CTRUNC.add_command_line_arguments('-mtzin', mtzin, '-mtzout', mtzout, '-colin', '/*/*/[{!s},{!s}]'.format(I,SIGI))
    # Run!
    CTRUNC.run()

    if not os.path.exists(mtzout):
        raise ExternalProgramError('CTRUNCATE has failed to convert intensities to SFs. {!s}\nOUT: {!s}\nERR: {!s}'.format(mtzin, CTRUNC.out, CTRUNC.err))

    return CTRUNC

def apply_rfree_set(refmtz, mtzin, mtzout):
    """Takes an input mtz and a reference mtz and transplants the rfree flags"""

    tempmtz = mtzin.replace('.mtz','.temp.mtz')

    refobj = MtzSummary(refmtz)
    newobj = MtzSummary(mtzin)

    # Copy data columns from new mtz file
    F1, SIGF1 = newobj.label.f, newobj.label.sigf
    if not (F1 and SIGF1):
        raise LabelError('No Amplitudes found in {!s}'.format(mtzin))

    # Get the spacegroup and the rfree from the reference mtz
    sgno = refobj.data.spacegroupno
    RFree2 = refobj.label.free
    if not sgno:
        raise LabelError('No Spacegroup information found in {!s}'.format(refmtz))
    if not RFree2:
        raise LabelError('No RFreeFlags found in {!s}'.format(refmtz))

    # Initialise Commander
    CAD = CommandManager('cad')
    # Set command arguments
    CAD.add_command_line_arguments('hklin1', mtzin, 'hklin2', refmtz, 'hklout', tempmtz)
    # Set inputs
    CAD.add_standard_input(['symmetry {!s}'.format(sgno),'labin file_number 1 E1={!s} E2={!s}'.format(F1, SIGF1), \
                                               'labout file_number 1 E1={!s} E2={!s}'.format(F1, SIGF1), \
                                               'labin file_number 2 E1={!s}'.format(RFree2), \
                                               'labout file_number 2 E1={!s}'.format(RFree2),'END'])
    # Run!
    CAD.run()

    if not os.path.exists(tempmtz):
        raise ExternalProgramError('CAD has failed to transplant RFree Flags. {!s}\nOUT: {!s}\nERR: {!s}'.format(mtzin, CAD.out, CAD.err))

    # Now complete the free r set

    # Initialise Commander
    FREE = CommandManager('freerflag')
    # Set command arguments
    FREE.add_command_line_arguments('hklin', tempmtz, 'hklout', mtzout)
    # Set inputs
    FREE.add_standard_input(['COMPLETE FREE={!s}'.format(RFree2), 'END'])
    # Run!
    FREE.run()

    if not os.path.exists(mtzout):
        raise ExternalProgramError('freerflag has failed to complete the RFree Flag set. {!s}\nOUT: {!s}\nERR: {s}'.format(tempmtz, CAD.out, CAD.err))

    os.remove(tempmtz)

    return CAD, FREE

