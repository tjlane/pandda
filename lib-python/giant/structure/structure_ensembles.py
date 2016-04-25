import os, sys, glob, time, re

class identical_structure_ensemble(object):
    """Class for collating and comparing multiple observations of the same structure"""

    def __init__(self, ref_hierarchy):
        """Initialise the comparison object"""

        self.data = multiple_data_collection()
        self.ref_hierarchy = ref_hierarchy

        # Set the collection ids to the residue labels
        residue_labels = [(rg.parent().id, rg.resid()) for rg in ref_hierarchy.residue_groups()]
        self.data.set_collection_ids(residue_labels)
        # List of just chain ids
        chain_ids = [rg.parent().id for rg in ref_hierarchy.residue_groups()]
        self.data.add_collection_info(info_name='chain ids', info_values=chain_ids)
        # List of just residue ids
        residue_ids = [rg.resid() for rg in ref_hierarchy.residue_groups()]
        self.data.add_collection_info(info_name='residue ids', info_values=residue_ids)
        # List of which residues are protein
        is_protein = [rg.parent().is_protein() for rg in ref_hierarchy.residue_groups()]
        self.data.add_collection_info(info_name='is protein', info_values=is_protein)
        # List of which residues are in multiple conformations
        has_conformers = [rg.have_conformers() for rg in ref_hierarchy.residue_groups()]
        self.data.add_collection_info(info_name='has conformers', info_values=has_conformers)
        # List of number of conformers per residue
        num_conformers = [len(rg.conformers()) for rg in ref_hierarchy.residue_groups()]
        self.data.add_collection_info(info_name='number of conformers', info_values=num_conformers)
        # List of number of unique resnames for each residue - XXX I'm not sure when this can be more than one?
        residue_names = [[s for s in rg.unique_resnames()] for rg in ref_hierarchy.residue_groups()]
        self.data.add_collection_info(info_name='residue names', info_values=residue_names)
        # List of residues types (amino, water, other, etc)
        residue_types = [iotbx.pdb.common_residue_names_get_class(rg.unique_resnames()[0]) for rg in ref_hierarchy.residue_groups()]
        self.data.add_collection_info(info_name='residue types', info_values=residue_types)
        # List of atom ids for each residue
        residue_atom_labels = [[a.id_str() for a in rg.atoms()] for rg in ref_hierarchy.residue_groups()]
        self.data.add_collection_info(info_name='atom labels', info_values=residue_atom_labels)

    def add_structures(self, new_hierarchies, hierarchy_ids=None, verbose=True):
        """Add hierarchies to the analysis object - Iterate through the residues and extract values from across the datasets"""

        # If no ids given, simply assign indexing values
        if not hierarchy_ids: hierarchy_ids = range(1, len(new_hierarchies)+1)

        # Report string to be returned
        report_string = []

        backbone_atom_names = [" CA ", " C  ", " O  ", " N  "]

        # Extract the reference residue types
        ref_residue_types = self.data.get_collection_info('residue types')
        ref_residue_names = self.data.get_collection_info('residue names')
        ref_chain_ids =     self.data.get_collection_info('chain ids')

        print('===================================>>>')
        print('Calculating B-factors for residues across the datasets')
        print('===================================>>>')

        # Atomic selection cache for the reference structure - should be the same across the structures
        selection_cache = self.ref_hierarchy.atom_selection_cache()

        # String of 1-letter codes for the protein
        amino_sequence = ''
        current_chain = ref_chain_ids[self.data.collection_ids[0]]

        for res_id in self.data.collection_ids:

            # Initialise a new object to hold the data for one residue
            res_collection = data_collection()
            res_collection.set_entry_ids(hierarchy_ids)

            # Boolean mask for the selected residue
            res_bool_selection = selection_cache.selection("chain '{!s}' and resid '{!s}'".format(*res_id))

            # Information to collect for the residues
            backbone_mean_bfactor = []
            sidechain_mean_bfactor = []
            total_mean_bfactor = []

            # Switches for tracking and reporting new chains
            if current_chain != ref_chain_ids[res_id]:
                report_string.append('\rAnalysing Chain {!s}: {!s}'.format(current_chain, amino_sequence))
                print(report_string[-1])
                report_string.append('\rChain {!s} Analysed.'.format(current_chain))
                print(report_string[-1])
                # Reset sequence and update chain
                amino_sequence = ''
                current_chain = ref_chain_ids[res_id]

            try:
                res_class = iotbx.pdb.common_residue_names_get_class(ref_residue_names[res_id][0])
                if res_class == 'common_amino_acid':
                    amino_sequence += iotbx.pdb.amino_acid_codes.one_letter_given_three_letter[ref_residue_names[res_id][0]]
                elif res_class == 'common_rna_dna':
                    amino_sequence += '<DNA/RNA>'
                elif res_class == 'ccp4_mon_lib_rna_dna':
                    amino_sequence += '<DNA/RNA>'
                elif res_class == 'common_water':
                    amino_sequence += 'o'
                elif res_class == 'common_small_molecule':
                    amino_sequence += '<MOL>'
                elif res_class == 'common_element':
                    amino_sequence += 'i'
                elif res_class == 'other':
                    amino_sequence += 'X'
                else:
                    raise Exception()
            except:
                amino_sequence += '?'

            if len(amino_sequence) >= 50:
                report_string.append('\rAnalysing Chain {!s}: {!s}...'.format(current_chain, amino_sequence))
                print(report_string[-1])
                amino_sequence = ''
            else:
                print('\rAnalysing Chain {!s}: {!s}-->'.format(current_chain, amino_sequence), end=''); sys.stdout.flush()

            # Iterate through the hierarchies and extract values for this residue
            for i_hierarchy in new_hierarchies:

                # Select the atoms for this residue
                new_root = i_hierarchy.select(res_bool_selection)

                # Calculate the sidechain and backbone b-factors for amino acids
                if ref_residue_types[res_id] == 'common_amino_acid':

                    # Pull out the backbone atoms
                    backbone_atoms = [at for at in new_root.atoms() if at.name in backbone_atom_names]
                    if not backbone_atoms: backbone_mean_bfactor.append(numpy.nan)
                    else:                  backbone_mean_bfactor.append(flex.mean(flex.double([at.b for at in backbone_atoms])))

                    # Pull out the sidechain atoms
                    sidechain_atoms = [at for at in new_root.atoms() if at.name not in backbone_atom_names]
                    if not sidechain_atoms: sidechain_mean_bfactor.append(numpy.nan)
                    else:                   sidechain_mean_bfactor.append(flex.mean(flex.double([at.b for at in sidechain_atoms])))

                else:
                    # If not an amino acid, just append None
                    backbone_mean_bfactor.append(numpy.nan)
                    sidechain_mean_bfactor.append(numpy.nan)

                # Calculate the mean b-factor for the whole residue
                total_atoms = [at for at in new_root.atoms()]
                if not total_atoms: total_mean_bfactor.append(numpy.nan)
                else:               total_mean_bfactor.append(flex.mean(flex.double([at.b for at in total_atoms])))

            res_collection.add_data(data_name='backbone_mean_bfactor', data_values=backbone_mean_bfactor)
            res_collection.add_data(data_name='sidechain_mean_bfactor', data_values=sidechain_mean_bfactor)
            res_collection.add_data(data_name='total_mean_bfactor', data_values=total_mean_bfactor)

            self.data.add_collection(collection_id=res_id, collection=res_collection)

        report_string.append('\rAnalysing Chain {!s}: {!s}   '.format(current_chain, amino_sequence))
        print(report_string[-1])
        report_string.append('\rChain {!s} Analysed.'.format(current_chain))
        print(report_string[-1])

        return '\n'.join(report_string)


