
import iotbx.pdb
from scitbx.math import dihedral_angle

from giant.structure.select import extract_backbone_atoms

########################################################################################

def get_all_phi_psi_for_hierarchy(pdb_h, deg=True):
    """Calculate all phi-psi angles for a structure"""
    for prev, current, next in generate_residue_triplets(pdb_h):
        phi, psi = calculate_residue_phi_psi_angles(prev=prev, current=current, next=next, deg=deg)
        yield (current, (phi,psi))

def generate_residue_triplets(pdb_h):
    """Generate groups of three residues (this is an iterator)"""
    for chain in pdb_h.chains():
        # Skip if not protein
        if not chain.is_protein(): continue
        for conf in chain.conformers():
            # Reset for each conformer
            prev = current = next = None
            for res in conf.residues():
                # Shuffle all variables down one
                prev = current
                current = next
                next = res
                # Do we have three residues?
                if not (prev and current and next):
                    continue
                # Are they continuous residues
                if not (prev.resseq_as_int()+1 == current.resseq_as_int() == next.resseq_as_int()-1):
                    continue
                # Return the residue triplet
                yield (prev, current, next)

########################################################################################

def calculate_residue_phi_psi_angles(prev, current, next, deg=True):
    """Calculate phi-psi angles for a set of three residues. Returns (phi, psi)."""
    phi = dihedral_angle(sites=extract_phi_sites(prev=prev, current=current), deg=deg)
    psi = dihedral_angle(sites=extract_psi_sites(current=current, next=next), deg=deg)
    return (phi, psi)

########################################################################################

def extract_psi_sites(current, next):
    """Extract sites to calculate psi angle from the selected residues"""
    return [a.xyz for a in extract_psi_atoms(current=current, next=next)]

def extract_phi_sites(prev, current):
    """Extract sites to calculate phi angle from the selected residues"""
    return [a.xyz for a in extract_phi_atoms(prev=prev, current=current)]

########################################################################################

def extract_psi_atoms(current, next):
    """Extract atoms to calculate phi angle from the selected residues"""
    n_curr = ca_curr = c_curr = n_next = None
    n_next, ca_next, c_next = extract_backbone_atoms(residue=next)
    n_curr, ca_curr, c_curr = extract_backbone_atoms(residue=current)
    return (n_curr, ca_curr, c_curr, n_next)

def extract_phi_atoms(prev, current):
    """Extract sites to calculate phi angle from the selected residues"""
    c_prev = n_curr = ca_curr = c_curr = None
    n_prev, ca_prev, c_prev = extract_backbone_atoms(residue=prev)
    n_curr, ca_curr, c_curr = extract_backbone_atoms(residue=current)
    return (c_prev, n_curr, ca_curr, c_curr)

########################################################################################

