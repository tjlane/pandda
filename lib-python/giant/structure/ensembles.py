from __future__ import print_function

import os, sys, glob, time, re

import numpy, pandas

import iotbx.pdb

from scitbx.array_family import flex

from bamboo.common import Meta, Info
from bamboo.plot.histogram import simple_histogram

from giant.structure import make_label
from giant.structure.dihedrals import get_all_phi_psi_for_hierarchy
from giant.structure.b_factors import BfactorStatistics, normalise_b_factors_to_z_scores
from giant.structure.select import non_h, backbone, sidechains

from giant.structure.iterators import residues_via_conformers, conformers_via_residue_groups

#def label_from_rg(rg):
#def label_from_res(res):

class StructureCollection(object):
    """Class to hold information about a set of related structures"""

    def __init__(self, hierarchy_inputs=None, hierarchies=None, inputs=None, labels=None):

        # Object to stores all of the structure objects
        self.structures = Info(['hierarchies','inputs','labels'])
        # Unpack and store structure objects
        if hierarchy_inputs:
            self.structures.hierarchies = [i.hierarchy for i in hierarchy_inputs]
            self.structures.inputs      = [i.input     for i in hierarchy_inputs]
        else:
            self.structures.hierarchies = hierarchies
            self.structures.inputs      = inputs
        # Unpack and store labels
        if not labels: labels = range(len(self.structures.hierarchies))
        self.structures.labels = labels

        #Initialise output tables
        self.initialise_tables()

    @classmethod
    def from_files(cls, filenames, labels=None):
        """Initialise the object from files"""
        if not labels: labels=filenames
        return cls(hierarchy_inputs=[iotbx.pdb.hierarchy.input(f) for f in filenames], labels=labels)

    def initialise_tables(self):
        """Populate the index of the tables"""

        # Object to contain all of the output data tables
        self.tables = Info(['structures','residues','atoms'])
        self._initialise_structure_table()
        self._initialise_residue_table()

    def _initialise_structure_table(self):
        """Initialise the tables.structures object with structure labels"""

        self.tables.structures = pandas.DataFrame( data = None,
                                                   index   = pandas.Index(data=self.structures.labels, name='structure label'),
                                                   columns = []      )

    def _initialise_residue_table(self):
        """Initialise the tables.residues object with residue labels"""

        residue_labels = [make_label(c) for c in conformers_via_residue_groups(self.structures.hierarchies[0])]
        self.tables.residues = pandas.Panel( data       = None,
                                             items      = pandas.Index(data=residue_labels, name=['chain','residue','altloc'], tupleize_cols=True),
                                             major_axis = pandas.Index(data=self.structures.labels, name='structure label'),
                                             minor_axis = pandas.Index(data=[], name='variables')     )

    def write_tables(self, out_dir):
        """Dump all data as csv files in out_dir"""

        if not os.path.exists(out_dir): os.mkdir(out_dir)

        # Write csv for each variable for residues
        for variable in self.tables.residues.minor_axis:
            filename = os.path.join(out_dir, variable+'.csv')
            print('Writing {}'.format(filename))
            self.tables.residues.xs(variable, axis=2).T.to_csv(filename)

        # Write csv for all global structure variables
        filename = os.path.join(out_dir, 'structures.csv')
        print('Writing {}'.format(filename))
        self.tables.structures.to_csv(filename)

    def write_graphs(self, out_dir):
        """Dump histograms in out_dir"""

        if not os.path.exists(out_dir): os.mkdir(out_dir)

        for variable in self.tables.residues.minor_axis:
            var_data = self.tables.residues.xs(variable, 2)
            for res_lab in self.tables.residues.items:
                filename = os.path.join(out_dir,'{}-{}.png'.format('-'.join(res_lab), variable))
                col_data = var_data[res_lab]
                if col_data.any() and (col_data.nunique()>1):
                    print('Writing {}'.format(filename))
                    simple_histogram( filename = filename,
                                      data = col_data,
                                      title = '{} histogram for {!s}'.format(variable, res_lab),
                                      x_lab = variable  )

        for variable in self.tables.structures.columns:
            col_data = self.tables.structures[variable]
            filename = os.path.join(out_dir,'{}.png'.format(variable))
            if col_data.any() and (col_data.nunique()>1):
                print('Writing {}'.format(filename))
                simple_histogram( filename = filename,
                                  data = col_data,
                                  title = 'histogram for {!s}'.format(variable),
                                  x_lab = variable  )

    def load_all(self):
        self.load_all_structure_data()
        self.load_all_residue_data()
    def load_all_structure_data(self):
        self.extract_structure_info()
        self.calculate_structure_mean_b_factors()
    def load_all_residue_data(self):
        self.extract_residue_info()
        self.calculate_residue_mean_b_factors()
        self.calculate_residue_mean_normalised_b_factors()
        self.calculate_residue_phi_psi_angles()

    def extract_structure_info(self):
        """Extract information from the input objects, if given"""

        if not self.structures.inputs: raise Exception('No input objects given')

        self.tables.structures['r-work'] = numpy.nan
        self.tables.structures['r-free'] = numpy.nan
        self.tables.structures['res-high'] = numpy.nan
        self.tables.structures['res-low']  = numpy.nan

        for lab_h, pdb_i in zip(self.structures.labels, self.structures.inputs):
            print('Extracting Global PDB Metrics: Structure: {}'.format(lab_h))
            r_info = pdb_i.get_r_rfree_sigma()
            self.tables.structures.set_value(lab_h, 'r-work', r_info.r_work)
            self.tables.structures.set_value(lab_h, 'r-free', r_info.r_free)
            self.tables.structures.set_value(lab_h, 'res-high', r_info.high)
            self.tables.structures.set_value(lab_h, 'res-low',  r_info.low)

    def calculate_structure_mean_b_factors(self):
        """Calculate global B-factors for a structure"""

        self.tables.structures['mean-b-all']       = numpy.nan
        self.tables.structures['mean-b-protein']   = numpy.nan
        self.tables.structures['mean-b-backbone']  = numpy.nan
        self.tables.structures['mean-b-sidechain'] = numpy.nan

        self.tables.structures['rms-b-all']       = numpy.nan
        self.tables.structures['rms-b-protein']   = numpy.nan
        self.tables.structures['rms-b-backbone']  = numpy.nan
        self.tables.structures['rms-b-sidechain'] = numpy.nan

        for lab_h, pdb_h in zip(self.structures.labels, self.structures.hierarchies):
            print('Calculating Global Mean B-Factors: Structure: {}'.format(lab_h))

            b_stats = BfactorStatistics.from_pdb(pdb_hierarchy=pdb_h.deep_copy())

            self.tables.structures.set_value(lab_h, 'mean-b-all',       b_stats.all.mean        )
            self.tables.structures.set_value(lab_h, 'mean-b-protein',   b_stats.protein.mean    )
            self.tables.structures.set_value(lab_h, 'mean-b-backbone',  b_stats.backbone.mean   )
            self.tables.structures.set_value(lab_h, 'mean-b-sidechain', b_stats.sidechain.mean  )

            self.tables.structures.set_value(lab_h, 'rms-b-all',       b_stats.all.biased_standard_deviation        )
            self.tables.structures.set_value(lab_h, 'rms-b-protein',   b_stats.protein.biased_standard_deviation    )
            self.tables.structures.set_value(lab_h, 'rms-b-backbone',  b_stats.backbone.biased_standard_deviation   )
            self.tables.structures.set_value(lab_h, 'rms-b-sidechain', b_stats.sidechain.biased_standard_deviation  )

    def extract_residue_info(self):
        """Extract information about each of the residues"""

        self.tables.residues.loc[:,:,'num_conformers'] = numpy.nan
        self.tables.residues.loc[:,:,'num_atoms']      = numpy.nan
        self.tables.residues.loc[:,:,'occupancy'] = numpy.nan

        for lab_h, pdb_h in zip(self.structures.labels, self.structures.hierarchies):
            print('Extracting Residue Info: Structure: {}'.format(lab_h))
            for c in conformers_via_residue_groups(pdb_h):
                res_lab = make_label(c)
                self.tables.residues.set_value(res_lab, lab_h, 'num_conformers', len(c.parent().conformers()))
                self.tables.residues.set_value(res_lab, lab_h, 'num_atoms', c.atoms_size())

                p_occ = [o for o in c.atoms().extract_occ() if o<1.0]
                if c.altloc and (len(set(p_occ)) == 1):
                    self.tables.residues.set_value(res_lab, lab_h, 'occupancy', p_occ[0])

    def calculate_residue_mean_b_factors(self):
        """Extract Mean-B values in each of the structures"""

        self.tables.residues.loc[:,:,'mean-b-all']       = numpy.nan
        self.tables.residues.loc[:,:,'mean-b-backbone']  = numpy.nan
        self.tables.residues.loc[:,:,'mean-b-sidechain'] = numpy.nan

        for lab_h, pdb_h in zip(self.structures.labels, self.structures.hierarchies):
            print('Calculating Local Mean B-Factors: Structure: {}'.format(lab_h))
            cache = pdb_h.atom_selection_cache()
            # Non-Hydrogens
            for c in conformers_via_residue_groups(non_h(hierarchy=pdb_h, cache=cache)):
                res_lab = make_label(c)
                res_mean_b = flex.mean_weighted(c.atoms().extract_b(), c.atoms().extract_occ())
                self.tables.residues.set_value(res_lab, lab_h, 'mean-b-all', res_mean_b)
            # Backbone Atoms
            for c in conformers_via_residue_groups(backbone(hierarchy=pdb_h, cache=cache)):
                res_lab = make_label(c)
                res_mean_b = flex.mean_weighted(c.atoms().extract_b(), c.atoms().extract_occ())
                self.tables.residues.set_value(res_lab, lab_h, 'mean-b-backbone', res_mean_b)
            # Sidechain Atoms
            for c in conformers_via_residue_groups(sidechains(hierarchy=pdb_h, cache=cache)):
                res_lab = make_label(c)
                res_mean_b = flex.mean_weighted(c.atoms().extract_b(), c.atoms().extract_occ())
                self.tables.residues.set_value(res_lab, lab_h, 'mean-b-sidechain', res_mean_b)

    def calculate_residue_mean_normalised_b_factors(self):
        """Extract Mean-B values in each of the structures"""

        self.tables.residues.loc[:,:,'mean-bz-all']       = numpy.nan
        self.tables.residues.loc[:,:,'mean-bz-backbone']  = numpy.nan
        self.tables.residues.loc[:,:,'mean-bz-sidechain'] = numpy.nan

        for lab_h, pdb_h in zip(self.structures.labels, self.structures.hierarchies):
            print('Calculating Local Normalised Mean B-Factors: Structure: {}'.format(lab_h))
            pdb_h_z = normalise_b_factors_to_z_scores(pdb_hierarchy=pdb_h, method='all')
            cache = pdb_h_z.atom_selection_cache()
            # Non-Hydrogens
            for c in conformers_via_residue_groups(non_h(hierarchy=pdb_h_z, cache=cache)):
                res_lab = make_label(c)
                res_mean_b = flex.mean_weighted(c.atoms().extract_b(), c.atoms().extract_occ())
                self.tables.residues.set_value(res_lab, lab_h, 'mean-bz-all', res_mean_b)
            # Backbone Atoms
            for c in conformers_via_residue_groups(backbone(hierarchy=pdb_h_z, cache=cache)):
                res_lab = make_label(c)
                res_mean_b = flex.mean_weighted(c.atoms().extract_b(), c.atoms().extract_occ())
                self.tables.residues.set_value(res_lab, lab_h, 'mean-bz-backbone', res_mean_b)
            # Sidechain Atoms
            for c in conformers_via_residue_groups(sidechains(hierarchy=pdb_h_z, cache=cache)):
                res_lab = make_label(c)
                res_mean_b = flex.mean_weighted(c.atoms().extract_b(), c.atoms().extract_occ())
                self.tables.residues.set_value(res_lab, lab_h, 'mean-bz-sidechain', res_mean_b)

    def calculate_normalised_normalised_residue_b_factors(self):
        """Calculate twice-normalised residue b-factors"""

        pass

    def calculate_residue_phi_psi_angles(self):
        """Extract phi-psi angles for each of the structures"""

        if not self.structures.hierarchies: raise Exception('No input objects given')

        # Add column labels
        self.tables.residues.loc[:,:,'phi'] = numpy.nan
        self.tables.residues.loc[:,:,'psi'] = numpy.nan

        for lab_h, pdb_h in zip(self.structures.labels, self.structures.hierarchies):
            print('Extracting Phi-Psi Angles: Structure: {}'.format(lab_h))
            # Note these are residue objects so model->chain->conformer->residue
            for res, (phi, psi) in get_all_phi_psi_for_hierarchy(pdb_h=pdb_h):
                # Create residue label - need to account for fact that all single confs are duplicated
                if res.is_pure_main_conf:
                    res_lab = (res.parent().parent().id, res.resid(), '')
                else:
                    res_lab = (res.parent().parent().id, res.resid(), res.parent().altloc)
                # Check to see if already done so only populate once
                if numpy.isnan(self.tables.residues.get_value(res_lab, lab_h, 'phi')):
                    self.tables.residues.set_value(res_lab, lab_h, 'phi', phi)
                # Populate values
                if numpy.isnan(self.tables.residues.get_value(res_lab, lab_h, 'psi')):
                    self.tables.residues.set_value(res_lab, lab_h, 'psi', psi)

    def write_normalised_b_factor_structures(self):
        """Write out the normalised b-factor structures"""

        for lab_h, pdb_h in zip(self.structures.labels, self.structures.hierarchies):
            print('Writing Normalised B-Factors: Structure: {}'.format(lab_h))
            pdb_h_z = normalise_b_factors_to_z_scores(pdb_hierarchy=pdb_h, method='backbone')
            pdb_h_z.write_pdb_file(file_name='{}-normalised-b.pdb'.format(lab_h))


