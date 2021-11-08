import iotbx.pdb

# TODO move to mulch.datasets TODO

def strip_pdb_to_string(input_pdb, remove_ter=False, remove_end=False):
    """Remove all ter cards from a pdb file. Output to string"""
    assert remove_ter or remove_end, 'No functions to perform'
    pdb_lines = open(str(input_pdb), 'r').readlines()
    out_lines = []
    for line in pdb_lines:
        if   remove_ter and (line[0:3]=='TER'):   pass
        elif remove_end and (line[0:3]=='END'):   pass
        else:   out_lines.append(line)
    return ''.join(out_lines)

def strip_pdb_to_input(input_pdb, remove_ter=False, remove_end=False):
    """Removes all of the TER and END records and regenerates them"""
    stripped = strip_pdb_to_string(input_pdb  = input_pdb,
                                        remove_ter = remove_ter,
                                        remove_end = remove_end     )
    return iotbx.pdb.hierarchy.input(pdb_string=stripped)

