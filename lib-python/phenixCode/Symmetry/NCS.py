#! /usr/local/python/python2.7.3-64bit/bin/python

print('IMPORTING FROM PHENIX!')

from Bamboo.Common.Command import CommandManager

def find_NCS_transforms(pdbfile, mtzfile=None):
    """usesi phenix.find_ncs to find the symmetry operators required to build the ASU"""

    PHENIX = CommandManager('phenix.find_ncs')
    PHENIX.SetArguments(pdbfile)
    PHENIX.Run()

    operators = []

    group_no = 0
    operator_no = 0

    output = PHENIX.out.split('\n')

    for i, line in enumerate(output):
        if line.startswith('GROUP'):
            new_group_no = int(line.lstrip('GROUP').strip())
            assert new_group_no == group_no + 1
            group_no = new_group_no

        if line.startswith('OPERATOR'):
            new_operator_no = int(line.lstrip('OPERATOR').strip())
            assert new_operator_no == operator_no + 1
            operator_no = new_operator_no

            rot1 = output[i+3]
            assert rot1.startswith('ROTA 1:')
            rot1 = rot1.lstrip('ROTA 1:').split()

            rot2 = output[i+4]
            assert rot2.startswith('ROTA 2:')
            rot2 = rot2.lstrip('ROTA 2:').split()

            rot3 = output[i+5]
            assert rot3.startswith('ROTA 3:')
            rot3 = rot3.lstrip('ROTA 3:').split()

            tran = output[i+6]
            assert tran.startswith('TRANS:')
            tran = tran.lstrip('TRANS:').split()

            op_matrix = rot1+rot2+rot3+tran

            operators.append((group_no, operator_no, op_matrix))

            print operators[-1]

    if not operators:
        print('NO NCS FOUND')

    return operators

