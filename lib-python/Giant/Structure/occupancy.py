from scitbx.array_family import flex

def normalise_occupancies(hierarchy, exclude_conformers=None, max_occ=1.0, min_occ=0.0, in_place=False):
    """Normalise the occupancies of a hierarchy so that the occupancies for a residue sum to 1.0"""

    if exclude_conformers is None: exclude_conformers = []

    # Edit original copy or create new one?
    if (not in_place): hierarchy = hierarchy.deep_copy()
    # Iterate through the output structure, and normalise the occupancies if necessary
    for rg in hierarchy.residue_groups():
        # No conformers - don't normalise
        if not rg.have_conformers(): continue

        # Get the groups to change and the groups to keep constant
        ag_chnge = [ag for ag in rg.atom_groups() if ag.altloc and (ag.altloc not in exclude_conformers)]
        if not ag_chnge: continue
        ag_const = [ag for ag in rg.atom_groups() if ag.altloc and (ag.altloc in exclude_conformers)]
        # Calculate the total occupancy of the groups
        occ_chnge = sum([max(ag.atoms().extract_occ()) for ag in ag_chnge])
        occ_const = sum([max(ag.atoms().extract_occ()) for ag in ag_const])

        if (occ_const+occ_chnge < min_occ):
            # Occupancy too low - normalise to minimum value
            new_occ_chnge = min_occ - occ_const
        elif (occ_const+occ_chnge > max_occ):
            # Occupancy too high - normalise to maximum value
            new_occ_chnge = max_occ - occ_const
        else:
            # Not selected for normalising
            continue

        # Normalise residue occupancies
        assert new_occ_chnge > 0.0, 'Occupancy of INVARIABLE AGs is more than 1.0: {}'.format(occ_const)
        [ag.atoms().set_occ(ag.atoms().extract_occ()*new_occ_chnge/occ_chnge) for ag in ag_chnge]

    return hierarchy

def set_conformer_occupancy(hierarchy, conf_id, conf_occ, in_place=False):
    """Normalise the occupancies of a hierarchy so that the occupancies for a residue sum to 1.0"""

    # Edit original copy or create new one?
    if (not in_place): hierarchy = hierarchy.deep_copy()
    # Iterate through the output structure, and set the occupancies where selected
    for rg in hierarchy.residue_groups():
        # Calculate the total occupancy for this residue group
        for ag in rg.atom_groups():
            if ag.altloc == conf_id: ag.atoms().set_occ(flex.double([conf_occ]*len(ag.atoms())))
    return hierarchy

