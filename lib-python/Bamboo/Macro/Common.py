#! /usr/local/python/python2.7.3-64bit/bin/python

from Bamboo.Macro.Molecule import MacroMol
from Bamboo.Utils.Rdkit.Smile import find_structure_matches
from Bamboo.Wrappers.PdbUtils import isolate_residue_by_res_id

def isolate_compound_from_file_by_smile(pdbin, pdbouttemplate, reference_smile, allow_part=True):
    """Isolates residues matching `reference_smile` from pdbin. Returns a list of files of the isolated residues based on pdbout"""

    if not pdbouttemplate.endswith('.pdb'):
        pdbouttemplate += '.pdb'
    pdbouttemplate = pdbouttemplate.replace('.pdb','.{!s}.pdb')

    # Load the molecule
    try:
        macro = MacroMol(pdbin)
    except MacroMolError:
        raise

    residues_to_remove = []

    # Get residues that match the input smile
    for res in macro.getResidues():
        if res.smile:
            full, part = find_structure_matches(res.smile, reference_smile)
            if full or (allow_part and part):
                residues_to_remove.append(res)

    if not residues_to_remove:
        raise Exception('Could Not Find Compound in File! ({!s} in {!s})'.format(reference_smile, pdbin))

    pdbout = []
    for i, res in enumerate(residues_to_remove):

        out = pdbouttemplate.format(i)
        pdbout.append(out)
        isolate_residue_by_res_id(pdbin, out, chain=res.chain, resnum=res.resnum, model='*', inscode=res.inscode)

    return pdbout

def extract_compound_smiles(pdbin, find_boring_compounds=False):
    """Takes a PDB file and returns smiles strings of all of the compounds in the file"""

    try:
        m = MacroMol(pdbin)
    except MacroMolError:
        raise

    cpds = m.getUnknowns()

    # ADD FILTERING STEP!
    if find_boring_compounds:
        cpds.extend(m.getSolvent())

    [c.getSmiles() for c in cpds]

    return [c.smile for c in cpds]
