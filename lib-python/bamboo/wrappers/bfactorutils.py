import os, shutil, tempfile

from bamboo.common.command import CommandManager

def extract_bfactor_statistics(pdbfile):
    """Analyse the b-factors in a pdb file and return a dictionary with all of the averages, rms deviations and Z-scores for the b-factors in a structure"""

    # Create a temporary file for the output summary table and pdbfile
    temp_handle, temp_path = tempfile.mkstemp(suffix='.table', prefix='baverage_')
    temp_handle2, temp_path2 = tempfile.mkstemp(suffix='.pdb', prefix='baverage_')

    BAVERAGE = CommandManager('baverage')
    BAVERAGE.SetArguments('XYZIN',pdbfile,'RMSTAB',temp_path,'XYZOUT',temp_path2)
    BAVERAGE.SetInput(['END'])
    BAVERAGE.Run()

    if not os.path.exists(temp_path):
        raise ExternalProgramError('BAVERAGE has failed to calculate b-factor summary statistics for {!s}'.format(pdbfile))

    # Process Table and Standard Out
    table_contents = open(temp_path, 'r').read().split('\n')

    if not table_contents:
        raise ExternalProgramError('BAVERAGE has failed to calculate b-factor summary statistics for {!s}'.format(pdbfile))
    else:
        os.remove(temp_path)
        os.remove(temp_path2)

    return BAVERAGE, table_contents
