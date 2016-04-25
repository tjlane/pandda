
from scitbx.array_family import flex

def extract_ligands_from_hierarchy(pdb_hierarchy):
    """Extracts ligands from hierarchy objects. Returns as new hierarchy objects"""

    # Filter out waters
    no_waters = pdb_hierarchy.select(pdb_hierarchy.atom_selection_cache().selection('not water'))
    # Get only hetatms
    hetero_only = no_waters.select(flex.bool([a.hetero for a in no_waters.atoms_with_labels()]))
    hetero_cache = hetero_only.atom_selection_cache()

    # Build chain-resnum combinations
    res_ids = []
    for ch in hetero_only.chains():
        for res in ch.residue_groups():
            res_ids.append((ch.id, res.resseq, res.icode))

    roots = []
    # Select each residue at a time
    for chn, seq, icode in res_ids:
        atm_sel = hetero_cache.selection("chain '{!s}' and resseq '{!s}' and icode '{!s}'".format(chn,seq,icode))
        new_root = hetero_only.select(atm_sel)
        roots.append(new_root)

    return roots
