from scitbx.array_family import flex
from mmtbx.command_line.super import extract_sequence_and_sites

def get_calpha_sites(input_obj, structure_obj):
    """Gets the coordinates of alpha carbons from the structure"""

    seq, sites, sites_bool = extract_sequence_and_sites(input_obj)

#    return flex.vec3_double([s for s,n in zip(structure_obj.sites_cart(),input_obj.atoms().extract_name()) if n.strip()=='C'])
    return sites


