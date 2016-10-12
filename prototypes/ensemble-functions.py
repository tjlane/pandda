#! pandda.python

import os, sys, glob

import matplotlib
matplotlib.interactive(0)
from matplotlib import pyplot

import numpy
import pandas
import iotbx.pdb
import scipy.linalg

from scitbx.array_family import flex

def calc_com():
    # Calculate centre of mass of the reference dataset
    pandda.log('===================================>>>', True)
    pandda.log('Calculating Reference Structure Centre of Mass', True)

    # Should be in the reference frame
    ref_structure = pandda.reference_dataset().hierarchy()
    backbone_atoms = ref_structure.select(ref_structure.atom_selection_cache().selection('pepnames and (name CA or name C or name O or name N)')).atoms()
    backbone_sites = backbone_atoms.extract_xyz()
    backbone_wghts = flex.double([1]*len(backbone_sites))

    centre_of_mass = backbone_sites.mean_weighted(weights=backbone_wghts)

def calculate_mean_structure_and_protein_masks(self, deviation_cutoff):
    """Calculate the average of all of the structures, and create masks for each protein where residues deviate from the mean by more than `deviation_cutoff`"""

    self.log('===================================>>>', True)
    self.log('Calculating Mean Structure', True)
    self.log('===================================>>>')

    # TODO Make this reject points until consensus

    # Pull all c-alpha sites for each structure
    all_sites = numpy.array([d.transform_to_reference(points=d.get_calpha_sites(), method='global') for d in self.datasets.mask(mask_name='rejected - total', invert=True)])
    # Calculate the mean x,y,z for each c-alpha
    mean_sites = numpy.mean(all_sites, axis=0)
    # Differences from the mean
    diff_sites = all_sites - mean_sites
    # Euclidean norms of the distances moved
    diff_norms = numpy.apply_along_axis(numpy.linalg.norm, axis=2, arr=diff_sites)

    # TODO MOVE THIS TO THE STRUCTURE VARIATION FUNCTION TODO

    # TODO CREATE A HIERARCHY FOR THE MEAN STRUCTURE (AND WITH MEAN NORMALISED B-FACTORS?)

    # Create a list of masks for large-moving c-alphas
    residue_deviation_masks = []
    # Iterate by dataset, masking if the deviation of the calpha in the dataset is more than `deviation_cutoff`
    for calpha_shifts in diff_norms:
        residue_deviation_masks.append([1 if shift > deviation_cutoff else 0 for shift in calpha_shifts])

    # Save the masks
    self._average_calpha_sites = flex.vec3_double(mean_sites)
    self._residue_deviation_masks = residue_deviation_masks

def analyse_ensemble(pdbs):

    # Load the ensemble into hierarchies and sort them
    hierarchies = [iotbx.pdb.hierarchy.input(p) for p in pdbs]
    hierarchies = sorted(hierarchies, key=lambda h: h.input.get_r_rfree_sigma().r_free)
    hierarchies = [h.hierarchy for h in hierarchies]

    # Select the lowest r_free as the reference
    reference = hierarchies[0]

    # Create pandas object to hold residue information (items=atom_labels, major_axis=structure_labels, minor_axis=variables)
    ats_panel = pandas.Panel(   data       = None,
                                items      = [at.pdb_label_columns() for at in reference.atoms()],
                                major_axis = range(len(hierarchies)),
                                minor_axis = ['B','X','Y','Z']          )

    # Populate AG panels
    for i_h, h in enumerate(hierarchies):
        print 'Processing Structure {}'.format(i_h)
        for i_ag, ag in enumerate(h.atom_groups()):
            # Extract atoms object
            ag_ats = ag.atoms()
            # Populate relevant columns
            ags_panel[ag.id_str()+ag.altloc]['All-Mean-B'][i_h] = flex.mean(ag_ats.extract_b())
            # Centre of Mass of Residue
            c_o_m = ag_ats.extract_xyz().mean()
            ags_panel[ag.id_str()+ag.altloc]['All-Mean-X'][i_h] = c_o_m[0]
            ags_panel[ag.id_str()+ag.altloc]['All-Mean-Y'][i_h] = c_o_m[1]
            ags_panel[ag.id_str()+ag.altloc]['All-Mean-Z'][i_h] = c_o_m[2]
            # Centre of Mass of Backbone
            # Centre of Mass of Sidechain
            # Dihedral Angles
            # Rotamer IDs?

    # Populate AT panels
    for i_h, h in enumerate(hierarchies):
        print 'Processing Structure {}'.format(i_h)
        for i_at, at in enumerate(h.atoms()):
            # Populate relevant columns
            ats_panel[at.pdb_label_columns()]['B'][i_h] = at.b
            ats_panel[at.pdb_label_columns()]['X'][i_h] = at.xyz[0]
            ats_panel[at.pdb_label_columns()]['Y'][i_h] = at.xyz[1]
            ats_panel[at.pdb_label_columns()]['Z'][i_h] = at.xyz[2]

    return hierarchies, ags_panel, ats_panel

def analyse_structure_variability_1(self):
    """Go through all of the datasets and collect lots of different structural characteristics of the datasets for identifying odd datasets"""

    self.log('===================================>>>', True)
    self.log('Populating Ensemble Summary', True)

    # Extract the hierarchy for each of the datasets
    hierarchies = [d.hierarchy() for d in self.datasets.mask(mask_name='rejected - total', invert=True)]
    hierarchy_ids = [d.tag for d in self.datasets.mask(mask_name='rejected - total', invert=True)]

    # Add structures to the ensemble summary (auto-processes them)
    ensemble_log = self.ensemble_summary.add_structures(new_hierarchies=hierarchies, hierarchy_ids=hierarchy_ids)
    # This will already have been printed to stdout, so hide it.
    self.log(ensemble_log, hide=True)


def analyse_alignment_rt_variability(pandda):
    """Look at all of the rotation matrices for the local alignments and calculate the rms between neighbours"""

    assert pandda.params.alignment.method in ['global', 'local']

    if pandda.params.alignment.method == 'global':
        pandda.log('GLOBAL ALIGNMENT SELECTED - NOT ANALYSING ROTATION MATRICES')
        return

    if pandda.args.output.plot_graphs:
        import matplotlib
        matplotlib.interactive(0)
        from matplotlib import pyplot

    rot_identity = scitbx.matrix.identity(3)

    # Select datasets to analyse
    used_datasets = pandda.datasets.mask(mask_name='rejected - total', invert=True)

    # Reference c_alpha labels
    ref_c_alpha_labels = sorted(used_datasets[0].local_alignment_transforms().keys())

    # Array to hold the output data
    num_datasets = len(used_datasets)
    num_pairs =  len(ref_c_alpha_labels)-1
    output_diffs = numpy.zeros((num_datasets, num_pairs, 2))

    # Iterate through the datasets and pull out the alignment matrices
    for d_num, d_handler in enumerate(used_datasets):

        alignments = d_handler.local_alignment_transforms()
        alignment_keys = sorted(alignments.keys())

        assert alignment_keys == ref_c_alpha_labels

        # Iterate through adjacent pairs of matrices
        for i in range(0, num_pairs):

            # Label and lsq fit for the current calpha
            calpha_1 = alignment_keys[i]
            assert calpha_1 == ref_c_alpha_labels[i]
            rt_1 = alignments[calpha_1].rt()
            # And for the next calpha
            calpha_2 = alignment_keys[i+1]
            assert calpha_2 == ref_c_alpha_labels[i+1]
            rt_2 = alignments[calpha_2].rt()

            # Calculate the mapping from one frame to the other
            rt_1_2 = rt_1 * rt_2.inverse()
            # Calculate the angle of the rotation matrix
            theta_rad = scitbx.math.math.acos((rt_1_2.r.trace()-1)/2.0)
            theta_deg = theta_rad * 180.0/scitbx.math.math.pi
            # Calculate the length of the shift
            t_shift =  rt_1_2.t.norm_sq()**0.5

            # Append to the array
            output_diffs[d_num, i, :] = theta_deg, t_shift

    # Directory to write the output to
    var_out_dir = pandda.output_handler.get_dir('analyses')

    # Write out graphs
    if pandda.args.output.plot_graphs:

        # Create labels
        labels = ['']*num_pairs
        for i in range(0, num_pairs, 5)+[num_pairs-1]:
            labels[i] = i+1
        # Clear the last n before the last one
        n = 4
        labels[-1-n:-1] = ['']*n

        # BOX PLOT OF ROTATION AND TRANSLATION SHIFTS
        fig = pyplot.figure()
        pyplot.title('ROTATION-TRANSLATION MATRIX VARIATION')
        # ADJACENT ANGLE VARIATION
        pyplot.subplot(2, 1, 1)
        pyplot.boxplot(x=output_diffs[:,:,0], notch=True, sym='.', widths=0.5, whis=[5,95], whiskerprops={'ls':'-'}, flierprops={'ms':1}, labels=labels) # whis='range'
        pyplot.xlabel('C-ALPHA')
        pyplot.ylabel('ANGLE CHANGE')
        # ADJACENT SHIFT VARIATION
        pyplot.subplot(2, 1, 2)
        pyplot.boxplot(x=output_diffs[:,:,1], notch=True, sym='.', widths=0.5, whis=[5,95], whiskerprops={'ls':'-'}, flierprops={'ms':1}, labels=labels) # whis='range'
        pyplot.xlabel('C-ALPHA')
        pyplot.ylabel('TRANSLATION CHANGE')
        # Apply tight layout to prevent overlaps
        pyplot.tight_layout()
        # Save both
        pyplot.savefig(os.path.join(var_out_dir, 'calpha_rt_variation.png'), format='png')
        pyplot.close(fig)

    # Write out to file
    numpy.savetxt(  fname = os.path.join(var_out_dir, 'calpha_rt_r_variation.csv'), X=output_diffs[:,:,0],
                    delimiter=',', newline='\n' )
    numpy.savetxt(  fname = os.path.join(var_out_dir, 'calpha_rt_t_variation.csv'), X=output_diffs[:,:,1],
                    delimiter=',', newline='\n' )

    # Do some other stuff...



if __name__ == '__main__':

    pdbs = glob.glob('*.pdb')
    out_dir = 'reduce_ensembles'

    from giant.structure.ensemble

    hierarchies, ags_panel, ats_panel = analyse_ensemble(pdbs=pdbs[0:10])

    # Calculate the PCAs of the coordinates of each atom
    for at_label in ats_panel.items:
        print at_label
        at_coords = ats_panel[at_label][['X','Y']]
        at_coords = at_coords.as_matrix()

        # 2 subplots sharing x-axis
        fig, (axis_1_1, axis_2_1) = pyplot.subplots(2, sharex=False)
        # 1st Plot - 1st Y-Axis
        line_1_1, = axis_1_1.plot(at_coords[:,0], at_coords[:,1])
        axis_1_1.set_xlabel('X')
        axis_1_1.set_ylabel('Y', color='k')
        # Joint legend
#        axis_1_1.legend(handles=[line_1_1, line_1_2, line_1_3], loc=4)
        # Title
        axis_1_1.set_title('Coord Plot')
        pyplot.savefig('./outfig1.png')
        pyplot.close(fig)

        # SVD
        U, s, Vt = scipy.linalg.svd(at_coords)
        V = Vt.T
#        vectors = (V[0], V[1], V[2])

        at_coords = at_coords.T

        # Manual
        cov_mat = numpy.cov(at_coords)
        eig_val_cov, eig_vec_cov = numpy.linalg.eig(cov_mat)


        print cov_mat


        break


