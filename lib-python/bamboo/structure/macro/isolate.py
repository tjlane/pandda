
from bamboo.rdkit_utils.smile import find_structure_matches
from bamboo.wrappers.pdbutils import isolate_residue_by_res_id

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
    for res in macro.get_residues():
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

