import unittest

import os,sys

import numpy, scipy.stats
from cctbx import uctbx, sgtbx
from scitbx.array_family import flex

from bamboo.common import Meta, Info
from giant.dataset import ElectronDensityMap
from pandda.analyse.classes import PanddaMapAnalyser, MapHolderList
from pandda.analyse import graphs

class TestMapAnalyser(unittest.TestCase):


    WRITE_OUTPUT = True

    def setUp(self):

        numpy.random.seed(100)

        # Define map setup
        self.map_size   = (10,10,10)
        self.map_max    = 10.0
        self.map_uc     = uctbx.unit_cell('10 10 10 90 90 90')
        self.map_sg     = sgtbx.space_group('P1')
        self.n_datasets = 50
        self.av_dataset_uncertainty = 0.1
        self.max_sadj_amp = 0.1
        self.sadj_beta  = 0.1

        # Sample random values for the test data
        self.true_map_vals = self.map_max*numpy.random.rand(*self.map_size)
        self.true_map_sadj = numpy.random.exponential(self.sadj_beta, size=self.map_size)/self.sadj_beta
        self.true_map_sadj = self.true_map_sadj * self.max_sadj_amp / max(self.true_map_sadj)

        # Create dataset uncertainty values
        self.dataset_uncs = numpy.random.binomial(1000, self.av_dataset_uncertainty, size=self.n_datasets)/1000.0

        #print 'True:', self.true_map_vals
        #print 'Sadj:', self.true_map_sadj
        #print 'Uncs:', self.dataset_uncs

        # Create values for the fake datasets
        self.dataset_maps = []
        for i in range(self.n_datasets):
            dataset_err = numpy.random.normal(scale=self.dataset_uncs[i], size=self.map_size)
            dataset_var = numpy.random.normal(scale=self.true_map_sadj)
            dataset_map = self.true_map_vals + dataset_err + dataset_var

            flex_map = flex.double(dataset_map.flatten())
            flex_map.reshape(flex.grid(self.map_size))

            for gp in flex.nested_loop(self.map_size):
                self.assertEqual(dataset_map[gp], flex_map[gp])

            dataset_obj = ElectronDensityMap(
                                map_data    = flex_map,
                                unit_cell   = self.map_uc,
                                map_indices = None,
                                map_size    = self.map_size,
                                sparse      = False     )
            dataset_obj.meta.num = i
            dataset_obj.meta.tag = str(i)
            self.dataset_maps.append(dataset_obj)

        self.map_list = MapHolderList()
        self.map_list.add(self.dataset_maps)

    def tear_down(self):

        pass

    def test_map_analyser(self):

        ma = PanddaMapAnalyser( dataset_maps     = self.map_list,
                                meta             = None,
                                statistical_maps = None,
                                parent           = None,
                                log              = None   )

        calc_mean_map = ma.calculate_mean_map(mask_name=None)
        calc_dset_unc = ma.calculate_map_uncertainties(masked_idxs=None, cpus=40)
        calc_maps_all = ma.calculate_statistical_maps(mask_name=None, cpus=40)
        calc_sadj_map = calc_maps_all.sadj_map
        calc_skew_map = calc_maps_all.skew_map
        calc_bimo_map = calc_maps_all.bimo_map

        print max(numpy.abs(calc_skew_map.data))
        print max(numpy.abs(calc_bimo_map.data))

        err_mean_map = numpy.abs(self.true_map_vals.flatten()   - calc_mean_map.data)
        err_sadj_map = numpy.abs(self.true_map_sadj.flatten()   - calc_sadj_map.data)
        err_dset_unc = numpy.abs(numpy.array(self.dataset_uncs) - numpy.array(calc_dset_unc))

        self.assertAlmostEqual(numpy.mean(calc_mean_map.data), self.map_max*0.5)

        self.assertGreater(scipy.stats.pearsonr(self.true_map_vals.flatten(),   calc_mean_map.data)[0],         0.98)
        self.assertGreater(scipy.stats.pearsonr(self.true_map_sadj.flatten(),   calc_sadj_map.data)[0],         0.98)
        self.assertGreater(scipy.stats.pearsonr(numpy.array(self.dataset_uncs), numpy.array(calc_dset_unc))[0], 0.95)

        self.assertAlmostEqual(max(err_mean_map), 0.107, places=3)
        self.assertAlmostEqual(max(err_dset_unc), 0.008, places=3)
        self.assertIs(calc_mean_map, ma.statistical_maps.mean_map)

        print numpy.array(calc_maps_all.sadj_map.data)

        if self.WRITE_OUTPUT:

            out_dir = '/work/pearce/test-map-analyser'
            if not os.path.exists(out_dir): os.mkdir(out_dir)

            map_pairs = [('mean_map.ccp4', calc_maps_all.mean_map),
                         ('sadj_map.ccp4', calc_maps_all.sadj_map),
                         ('skew_map.ccp4', calc_maps_all.skew_map),
                         ('bimo_map.ccp4', calc_maps_all.bimo_map)]

            for mn, mo in map_pairs:
                mo.to_file(filename=os.path.join(out_dir, mn), space_group=self.map_sg)
