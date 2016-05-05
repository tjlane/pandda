import os, shutil, tempfile

from bamboo.common.command import CommandManager

def get_pdb_bfactor_summary(pdb_file):
    """Analyse the b-factors in a pdb file and return a dictionary with all of the averages, rms deviations and Z-scores for the b-factors in a structure"""

    program, contents = extract_bfactor_statistics(pdb_file)
    # Process the output table from baverage
    processed_contents = process_baverage_output(contents)
    # Process the stdout from the program
    processed_stdout = process_baverage_stdout(program.out)

    return program, processed_contents

def process_baverage_output(contents):
    """Process baverage output"""

    chain = ''
    resid = ''

    processed = {}

    for line in contents:
        if not line:
            continue
        elif line.strip().startswith('#'):
            chain = line.split()[2]
            assert len(chain)==1, 'Chain ID is more than one letter!'
        else:
            if not chain:
                raise Exception('No Chain ID has been assigned!')

            resid, m_av, m_rms, s_av, s_rms, a_av, a_rms = map(float, line.split())
            resid = int(resid)

            processed[(chain, resid)] = bFactorSummary(chain, resid, m_av, m_rms, s_av, s_rms, a_av, a_rms)

    return processed

def extract_bfactor_statistics(pdb_file):
    """Analyse the b-factors in a pdb file and return a dictionary with all of the averages, rms deviations and Z-scores for the b-factors in a structure"""

    # Create a temporary file for the output summary table and pdb_file
    temp_handle, temp_path = tempfile.mkstemp(suffix='.table', prefix='baverage_')
    temp_handle2, temp_path2 = tempfile.mkstemp(suffix='.pdb', prefix='baverage_')

    BAVERAGE = CommandManager('baverage')
    BAVERAGE.add_command_line_arguments('XYZIN',pdb_file,'RMSTAB',temp_path,'XYZOUT',temp_path2)
    BAVERAGE.add_standard_input(['END'])
    BAVERAGE.run()

    if not os.path.exists(temp_path):
        raise ExternalProgramError('BAVERAGE has failed to calculate b-factor summary statistics for {!s}'.format(pdb_file))

    # Process Table and Standard Out
    table_contents = open(temp_path, 'r').read().split('\n')

    if not table_contents:
        raise ExternalProgramError('BAVERAGE has failed to calculate b-factor summary statistics for {!s}'.format(pdb_file))
    else:
        os.remove(temp_path)
        os.remove(temp_path2)

    return BAVERAGE, table_contents
