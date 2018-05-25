import os, shutil

import pytest

import math, numpy
from libtbx.utils import Failure, Sorry

from pandda.analyse import pandda_phil, run
from pandda.resources.test_data import TEST_DATA

@pytest.fixture(scope="session")
def short_pandda():

    orig_dir = os.getcwd()
    work_dir = TEST_DATA.extract_to_temporary_directory(choice=TEST_DATA.keys.BAZ2B_TEST_DATA)
    os.chdir(work_dir)
    data_dir = './data'
    assert os.path.exists(data_dir)

    import multiprocessing
    cpus = max(1, multiprocessing.cpu_count()-1)

    # Input options
    options_phil = """
    data_dirs='./data/*'
    pdb_style='*.dimple.pdb'

    cpus={cpus}

    high_res_upper_limit=2.0
    high_res_lower_limit=2.0
    dynamic_res_limits=False
    min_build_datasets=1

    ignore_datasets='BAZ2BA-x432'
    exclude_from_characterisation='BAZ2BA-x434'
    exclude_from_z_map_analysis='BAZ2BA-x430,BAZ2BA-x431'
    ground_state='BAZ2BA-x430,BAZ2BA-x433,BAZ2BA-x435,BAZ2BA-x436,BAZ2BA-x437,BAZ2BA-x439'

    """.format(cpus=cpus)

    # Parse input args
    cmd_interpr = pandda_phil.command_line_argument_interpreter()
    working_phil = pandda_phil.fetch(sources=[cmd_interpr.process(options_phil)])

    # Make pandda object and run setup
    yield run(params=working_phil.extract())

    shutil.rmtree(work_dir)
    os.chdir(orig_dir)

class TestShortPanDDA(object):

    def test_grid(self, short_pandda):
        grid = short_pandda.grid

        assert grid.grid_size()     == (136, 134, 116)
        assert grid.grid_size_1d()  == 2113984

        assert grid.cart_origin()   == pytest.approx((-3.753, -1.563, -0.825))
        assert grid.cart_size()     == pytest.approx((68.0, 67.0, 58.0))

        assert grid.global_mask().total_size()   == 289603
        assert grid.symmetry_mask().total_size() == 519709
        assert grid.global_mask().extent()       == ((7, 88, 46), (128, 37, 34))
        assert grid.symmetry_mask().extent()     == ((0, 0, 8), (135, 133, 90))

        for v1, v2 in [((5,10,15),(17.506, 23.126, 31.65))]:
            assert grid.cart2grid([v1])[0]  == pytest.approx(v2)
            assert grid.grid2cart([v2])[0]  == pytest.approx(v1)

        assert grid.grid_point_volume()                      == pytest.approx(grid.grid_spacing()**3)
        assert grid.unit_cell().parameters()                 == pytest.approx((68.0, 67.0, 58.0, 90.0, 90.0, 90.0))
        assert grid.space_group().info().symbol_and_number() == 'P 1 (No. 1)'

        for i in [1,5,9,34,78]:
            p = grid.partition.sites_cart[i]
            assert grid.partition.query_by_cart_points([p]) == [i]
            assert grid.partition.query_by_grid_points([map(int,v) for v in grid.cart2grid([p])]) == [i]

    def test_dataset_info(self, short_pandda):
        table = short_pandda.tables.dataset_info

        assert numpy.unique(table['space_group'].dropna().values).tolist() == ['C 2 2 21']
        for col, vals in [('uc_alpha',                          [  90.0,    90.0,    90.0,    90.0,    90.0,    90.0,    90.0,    90.0,    90.0  ]),
                          ('uc_beta',                           [  90.0,    90.0,    90.0,    90.0,    90.0,    90.0,    90.0,    90.0,    90.0  ]),
                          ('uc_gamma',                          [  90.0,    90.0,    90.0,    90.0,    90.0,    90.0,    90.0,    90.0,    90.0  ]),
                          ('uc_vol',                            [ 460589.0,   458423.0,  459145.0,  455473.0,  463276.0,  461747.0,  457984.0,  459525.0,  460928.0 ]),
                          ('uc_a',                              [  82.128,  82.101,  82.03,   81.403,  81.879,  81.851,  81.676,  82.307,  82.323]),
                          ('uc_b',                              [  96.768,  96.536,  96.69,   96.831,  97.194,  97.145,  96.902,  96.516,  96.598]),
                          ('uc_c',                              [  57.955,  57.84,   57.889,  57.784,  58.214,  58.071,  57.866,  57.846,  57.962]),
                          ('applied_b_factor_scaling',          [  -0.001,  -0.293,  -0.197,  -0.0,    -0.649,  -1.503,  -0.785,   0.227,   0.286]),
                          ('high_resolution',                   [   1.786,   1.772,   1.67,    1.655,   1.858,   1.946,   1.877,   1.933,   1.805]),
                          ('low_resolution',                    [  62.616,  62.541,  62.552,  62.31,   62.62,   62.595,  62.451,  62.627,  62.656]),
                          ('r_free',                            [   0.21,    0.211,   0.209,   0.204,   0.205,   0.21,    0.211,   0.232,   0.202]),
                          ('r_work',                            [   0.182,   0.185,   0.184,   0.178,   0.176,   0.172,   0.178,   0.191,   0.179]),
                          ('rmsd_to_reference',                 [   0.086,   0.076,   0.092,   0.0,     0.062,   0.082,   0.083,   0.091,   0.066]),
                          ('scaled_wilson_B',                   [  10.322,  10.375,  10.398,  10.663,  10.853,  10.626,  10.180,  10.810,  10.298]),
                          ('scaled_wilson_ln_dev',              [  89.948,  96.668,  67.625,   0.0,   179.468, 108.934,  91.579, 315.795,  84.694]),
                          ('scaled_wilson_ln_dev_z',            [  -0.063,   0.198,  -0.931,  -3.559,   3.416,   0.674,   0.0,     8.714,  -0.268]),
                          ('scaled_wilson_ln_rmsd',             [   0.069,   0.069,   0.045,   0.0,     0.163,   0.128,   0.066,   0.256,   0.057]),
                          ('scaled_wilson_ln_rmsd_z',           [   0.0,     0.0,    -0.675,  -1.939,   2.642,   1.658,  -0.084,   5.255,  -0.337 ]),
                          ('scaled_wilson_rmsd_<4A',            [  3980.166,  6016.336,  2924.713, 0.0,  1297.052,  1877.977,  2249.678,  6119.786,  5245.620]),
                          ('scaled_wilson_rmsd_>4A',            [ 24045.415, 16314.626, 18710.202, 0.0,  1236.838,  3629.565,  6116.451, 29972.906, 24577.607]),
                          ('scaled_wilson_rmsd_all',            [  4015.332,  6025.885,  2953.761, 0.0,  1297.022,  1879.261,  2253.270,  6154.856,  5273.026]),
                          ('scaled_wilson_rmsd_<4A_z',          [   0.437,   1.281,   0.0,    -1.212,  -0.675,  -0.434,  -0.280,   1.324,   0.962]),
                          ('scaled_wilson_rmsd_>4A_z',          [   0.511,   0.0,     0.158,  -1.079,  -0.997,  -0.839,  -0.675,  0.903,   0.547]),
                          ('scaled_wilson_rmsd_all_z',          [   0.432,   1.251,   0.0,    -1.203,  -0.675,  -0.437,  -0.285,   1.303,   0.944]),
                          ('unscaled_wilson_B',                 [  10.323,  10.667,  10.596,  10.663,  11.502,  12.129,  10.965,  10.583,  10.012]),
                          ('unscaled_wilson_ln_dev',            [  90.156, 267.592, 187.671,   0.0,   627.811,1071.687, 606.765, 306.806, 192.784]),
                          ('unscaled_wilson_ln_dev_z',          [  -0.675,   0.0,    -0.304,  -1.017,   1.369,   3.057,   1.289,   0.149,  -0.284]),
                          ('unscaled_wilson_ln_rmsd',           [   0.07,    0.164,   0.107,   0.0,     0.377,   0.596,   0.334,   0.265,   0.102]),
                          ('unscaled_wilson_ln_rmsd_z',         [  -0.628,   0.0,    -0.381,  -1.095,   1.422,   2.885,   1.135,   0.675, -0.414]),
                          ('unscaled_wilson_rmsd_<4A',          [  3977.052,  4508.280,  3651.870, 0.0,  9996.903,  23069.640, 11878.436,  7934.614,  9076.754]),
                          ('unscaled_wilson_rmsd_>4A',          [ 24091.098, 26670.960, 25673.839, 0.0, 24854.683,  55153.129, 34010.256, 22077.427, 14402.213]),
                          ('unscaled_wilson_rmsd_all',          [  4012.384,  4546.437,  3695.816, 0.0, 10009.845,  23096.820, 11899.792,  7947.976,  9080.197]),
                          ('unscaled_wilson_rmsd_<4A_z',        [  -0.677,  -0.586,  -0.732,  -1.357,   0.353,   2.588,   0.675,   0.0,     0.195]),
                          ('unscaled_wilson_rmsd_>4A_z',        [  -0.185,   0.441,   0.199,  -6.036,   0.0,     7.358,   2.224,  -0.675,  -2.539]),
                          ('unscaled_wilson_rmsd_all_z',        [  -0.675,  -0.583,  -0.729,  -1.362,   0.353,   2.596,   0.677,   0.0,     0.194]),
                         ]:
            for v1, v2 in zip(table[col].dropna().values.tolist(), vals):
                assert pytest.approx(v1, abs=1e-3) == v2

    def test_dataset_labelling(self, short_pandda):
        table = short_pandda.datasets.all_masks().table

        # Check correct datasets are loaded
        all_d = 'BAZ2BA-x430,BAZ2BA-x431,BAZ2BA-x433,BAZ2BA-x434,BAZ2BA-x435,BAZ2BA-x436,BAZ2BA-x437,BAZ2BA-x438,BAZ2BA-x439'
        assert table.index.values.tolist() == all_d.split(',')
        # Check ignored dataset was not loaded
        ign_d = 'BAZ2BA-x432'
        assert ign_d not in table.index

        # Check labels are as expected
        for col, sel in [('rejected - total',               ''),
                         ('noisy zmap',                     ''),
                         ('analysed',                       'BAZ2BA-x433,BAZ2BA-x434,BAZ2BA-x435,BAZ2BA-x436,BAZ2BA-x437,BAZ2BA-x438,BAZ2BA-x439'),
                         ('interesting',                    'BAZ2BA-x434'),
                         ('exclude_from_z_map_analysis',    'BAZ2BA-x430,BAZ2BA-x431'),
                         ('exclude_from_characterisation',  'BAZ2BA-x431,BAZ2BA-x434,BAZ2BA-x438'),
                         ('old datasets',                   ''),
                         ('new datasets',                   all_d),
                         ('reprocess',                      ''),
                         ('valid - all',                    all_d),
                         ('valid - old',                    ''),
                         ('valid - new',                    all_d),
                         ('characterisation',               'BAZ2BA-x430,BAZ2BA-x433,BAZ2BA-x435,BAZ2BA-x436,BAZ2BA-x437,BAZ2BA-x439'),
                         ('analyse - old',                  ''),
                         ('analyse - new',                  'BAZ2BA-x433,BAZ2BA-x434,BAZ2BA-x435,BAZ2BA-x436,BAZ2BA-x437,BAZ2BA-x438,BAZ2BA-x439'),
                         ('analyse - all',                  'BAZ2BA-x433,BAZ2BA-x434,BAZ2BA-x435,BAZ2BA-x436,BAZ2BA-x437,BAZ2BA-x438,BAZ2BA-x439'),
                         ('not for analysis',               'BAZ2BA-x430,BAZ2BA-x431'),
                         ('characterisation @ 2.0A',        'BAZ2BA-x430,BAZ2BA-x433,BAZ2BA-x435,BAZ2BA-x436,BAZ2BA-x437,BAZ2BA-x439'),
                         ('analysis @ 2.0A',                'BAZ2BA-x433,BAZ2BA-x434,BAZ2BA-x435,BAZ2BA-x436,BAZ2BA-x437,BAZ2BA-x438,BAZ2BA-x439'),
                         ('loading @ 2.0A',                 'BAZ2BA-x430,BAZ2BA-x433,BAZ2BA-x434,BAZ2BA-x435,BAZ2BA-x436,BAZ2BA-x437,BAZ2BA-x438,BAZ2BA-x439'),
                        ]:
            sel = table.index.isin(sel.split(','))
            assert table.loc[sel][col].all()
            assert not table.loc[~sel][col].any()

    def test_dataset_map_info(self, short_pandda):
        table = short_pandda.tables.dataset_map_info

        nan_cols = ['analysed_resolution','map_uncertainty','obs_map_mean','obs_map_rms','scl_map_mean','scl_map_rms','z_map_kurt','z_map_mean','z_map_skew','z_map_std']
        assert table.loc['BAZ2BA-x430'][nan_cols].apply(math.isnan).values.all()
        assert table.loc['BAZ2BA-x431'].apply(math.isnan).values.all()

        for col, vals in [('analysed_resolution',  [ 2.0,    2.0,    2.0,    2.0,    2.0,    2.0,    2.0  ] ),
                          ('map_uncertainty',      [ 0.128,  0.243,  0.164,  0.138,  0.109,  0.226,  0.147] ),
                          ('obs_map_mean',         [ 0.051,  0.047,  0.044,  0.044,  0.05 ,  0.061,  0.059] ),
                          ('obs_map_rms',          [ 1.2  ,  1.194,  1.199,  1.199,  1.197,  1.192,  1.197] ),
                          ('z_map_kurt',           [ 2.263,  3.638,  2.428,  2.414,  2.465,  3.109,  2.575] ),
                          ('z_map_mean',           [-0.036, -0.   , -0.008, -0.006,  0.   , -0.026, -0.017] ),
                          ('z_map_skew',           [ 0.027,  0.091,  0.002, -0.003,  0.001,  0.043,  0.048] ),
                          ('z_map_std',            [ 0.825,  0.987,  0.871,  0.801,  0.794,  0.905,  0.848] ),
                          ('map_uncertainty-2.0A', [ 0.123,  0.128,  0.243,  0.164,  0.138,  0.109,  0.226,  0.147] ),
                          ('obs_map_mean-2.0A',    [ 0.041,  0.051,  0.047,  0.044,  0.044,  0.05 ,  0.061,  0.059] ),
                          ('obs_map_rms-2.0A',     [ 1.198,  1.2  ,  1.194,  1.199,  1.199,  1.197,  1.192,  1.197] ),
                          ('scl_map_mean-2.0A',    [ 0.057,  0.043,  0.047,  0.045,  0.046,  0.048,  0.041,  0.045] ),
                          ('scl_map_rms-2.0A',     [ 1.172,  1.177,  1.194,  1.185,  1.18 ,  1.181,  1.177,  1.179] ),
                          ('scl_map_mean',         [ 0.043,  0.047,  0.045,  0.046,  0.048,  0.041,  0.045] ),
                          ('scl_map_rms',          [ 1.177,  1.194,  1.185,  1.18 ,  1.181,  1.177,  1.179] ),
                          ]:
            assert table[col].dropna().values.tolist() == pytest.approx(vals, abs=1e-3)

    def test_events(self, short_pandda):
        table = short_pandda.tables.event_info

        for col, vals in [('1-BDC',             [ 0.19 ]),
                          ('cluster_size',      [ 160 ]),
                          ('site_idx',          [ 1 ]),
                          ('x',                 [ 13.247 ]),
                          ('y',                 [ 42.437 ]),
                          ('z',                 [ 31.675 ]),
                          ('z_mean',            [ 3.66 ]),
                          ('z_peak',            [ 5.96 ]),
                          ('global_correlation_to_average_map', [ 0.7788 ]),
                          ('local_correlation_to_average_map',  [ 0.3571 ]),
                         ]:
            assert table[col].values.tolist() == pytest.approx(vals, abs=1e-3)

