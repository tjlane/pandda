import numpy
from scitbx.array_family import flex

def pairwise_dists(xyz1, xyz2):
    """Calculate pairwise distances between xyz1 and xyz2"""
    return numpy.array([flex.sqrt((xyz2 - c1).dot()) for c1 in xyz1])

