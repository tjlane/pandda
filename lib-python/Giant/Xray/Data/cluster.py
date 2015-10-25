
import itertools
import numpy, pandas

import iotbx.mtz
import libtbx.cluster
import scipy.cluster.hierarchy

from Giant.Xray.Data import crystalSummary, unitCellVariation
from Giant.Xray.unit_cell import lcv_from_unit_cells, pairwise_lcv, pairwise_ecv

class crystalGroup(object):
    def __init__(self, crystals):
        """Class to hold multiple related crystals"""
        self.crystals = crystals
        self.uc_stats = unitCellVariation([c.unit_cell for c in self.crystals])
        self.space_groups = list(set([c.space_group.info().symbol_and_number() for c in self.crystals]))
    def size(self):
        return len(self.crystals)
    def cluster(self):
        cluster_func = lambda a,b: lcv_from_unit_cells(a.unit_cell, b.unit_cell)
        return libtbx.cluster.HierarchicalClustering(self.crystals, cluster_func)
    def dendrogram(self, fname, link_func=pairwise_lcv):
        """Print a dendrogram of the variation in the unit cells"""
        return unit_cell_dendrogram(fname=fname,
                                    unit_cells=[c.unit_cell for c in self.crystals],
                                    link_func=link_func,
                                    labels=[c.id for c in self.crystals])
    ######################################################
    @classmethod
    def by_space_group(cls, crystals):
        sg_func = lambda c: c.space_group.info().symbol_and_number()
        return [cls(list(x[1])) for x in itertools.groupby(sorted(crystals, key=sg_func), key=sg_func)]
    @classmethod
    def by_unit_cell(cls, crystals):
        pass

def dendrogram(fname, link_mat, labels=None):
    from matplotlib import pyplot
    fig = pyplot.figure()
    ax1 = pyplot.subplot(1,1,1)
    dend = scipy.cluster.hierarchy.dendrogram(link_mat)
    # Change labels if requested
    if labels: ax1.set_xticklabels([labels[i] for i in dend['leaves']])
    # Make sure the labels are rotated
    xlocs, xlabels = pyplot.xticks()
    pyplot.setp(xlabels, rotation=90)
    for i, d in zip(dend['icoord'], dend['dcoord']):
        x = 0.5 * sum(i[1:3])
        y = d[1]
        if y < 0.5: continue
        pyplot.plot(x, y, 'ro')
        pyplot.annotate("%.3g" % y, (x, y), xytext=(0, -8), textcoords='offset points', va='top', ha='center')
    pyplot.tight_layout()
    fig.savefig(fname)
    return fig

def unit_cell_dendrogram(fname, unit_cells, link_func, labels=None):
    """Calculate a dendrogram for a set of mtz files/object, clustering by unit_cell parameters"""
    dist_mat = link_func(unit_cells=unit_cells)
    link_mat = scipy.cluster.hierarchy.linkage(dist_mat)
    dendrogram(fname=fname, link_mat=link_mat, labels=labels)

def sort_mtzs_by_spacegroup(mtz_files=None, mtz_objects=None, labels=None):
    """Group the mtzfiles by spacegroup and then cluster them by unit cell variation"""

    assert [mtz_files, mtz_objects].count(None) == 1, 'Provide mtz_files OR mtz_objects'
    if mtz_files: mtz_objects = [None]*len(mtz_files)
    else:         mtz_files   = [None]*len(mtz_objects)

    if labels: assert len(labels) == len(mtz_files), 'labels must be the same length as mtz_files/mtz_objects'
    else:      labels = range(len(mtz_files))

    # Create summary objects
    mtz_summaries = [crystalSummary(mtz_file=m_f, mtz_object=m_o, id=lab) for m_f,m_o,lab in zip(mtz_files, mtz_objects, labels)]
    # Group by SpaceGroup
    cgs = crystalGroup.by_space_group(mtz_summaries)

    return cgs

def cluster_crystals(crystal_group):

    pass

