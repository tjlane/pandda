
import os, glob

from scipy.stats import zscore

import numpy

from scitbx.array_family import flex
from scitbx.math.distributions import normal_distribution

from PANDDAs.Main import dataset_handler

pdb_file = '/home/npearce/Analyses/FragmentSoaks/BAZ2BA-Feb13-Soak/reference.pdb'
mtz_file = '/home/npearce/Analyses/FragmentSoaks/BAZ2BA-Feb13-Soak/reference.mtz'

dirs = sorted(glob.glob('/home/npearce/Analyses/FragmentSoaks/BAZ2BA-Feb13-Soak/ProcessedFragmentSoak/*/1-apo'))

map_sigmas = []

for i, ddir in enumerate(dirs):

    try:
        pdb_file = glob.glob(os.path.join(ddir, 'apo-BAZ2BA-*-refmac.pdb'))[0]
        mtz_file = glob.glob(os.path.join(ddir, 'apo-BAZ2BA-*-refmac.mtz'))[0]
    except:
        continue

#    print pdb_file
#    print mtz_file
#    print '==============>'

    dh = dataset_handler(-1, pdb_filename=pdb_file, mtz_filename=mtz_file)

#    # OBSERVED MAPS
#    obs_map = dh.create_fft_map(map_type='2mFo-DFc', d_min=None)
#    obs_map.apply_volume_scaling()
#    obs_mh = dh.create_map_handler()
#    obs_vals = obs_map.real_map_unpadded(False).as_1d()
#    sorted_obs_vals = flex.double(sorted(obs_vals))

    # DIFFERENCE MAPS
    diff_map = dh.create_fft_map(map_type='mFo-DFc', d_min=None)
    diff_map.apply_volume_scaling()
    diff_mh = dh.create_map_handler()
    diff_vals = diff_map.real_map_unpadded(False).as_1d()
    sorted_diff_vals = flex.double(sorted(diff_vals))

#    print 'OBS MAP: '
#    obs_map.statistics().show_summary()
#    print '==============>'
#    print 'DIFF MAP: '
#    diff_map.statistics().show_summary()
#    print '==============>'

    # Normal distribution and expected quantiles
    normal_dist = normal_distribution()
    all_qs = normal_dist.quantiles(len(sorted_diff_vals))

    # Get the quantile indices between -1.5, 1.5
    low  = (all_qs < 1.5).iselection()
    high = (all_qs > -1.5).iselection()
    midd = low.intersection(high)

#    print 'TOTAL: ', len(sorted_obs_vals)
#    print '  LOW: ', len(low),  len(sorted_obs_vals) - len(low)
#    print ' HIGH: ', len(high), len(sorted_obs_vals) - len(high)
#    print ' COMB: ', len(midd), len(sorted_obs_vals) - len(midd)
#    print '==============>'

    # Select the values from the sorted list for the expected quantiles
    midd_sorted_vals = sorted_diff_vals.select(midd)
    midd_quantiles = all_qs.select(midd)

    # Calculate the slope of the line from the expected quantiles
    params = numpy.polyfit(x=midd_quantiles, y=midd_sorted_vals, deg=1)

#    print 'TOTAL SIGMA: ', diff_map.statistics().sigma()
#    print '  EST SIGMA: ', params[0]

    print '{!s}%'.format(int(100.0*i/len(dirs))),'-',pdb_file[84:88],'-',round(params[0],3)
    print 'TOTAL SIGMA: ', round(diff_map.statistics().sigma(),3)

    map_sigmas.append((pdb_file[84:88],params[0]))

################################

print 'MEAN SIGMA: ', numpy.mean(zip(*map_sigmas)[1])
print ' STD SIGMA: ', numpy.std(zip(*map_sigmas)[1])

zs = zscore(zip(*map_sigmas)[1])
