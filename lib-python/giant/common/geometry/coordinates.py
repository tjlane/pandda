import numpy as np

def distance_between(point_1, point_2):

    p1 = np.array(point_1)
    p2 = np.array(point_2)

    assert p1.ndim == 1
    assert p1.ndim == 1

    diff = (p1 - p2)

    dist = np.sqrt((diff*diff).sum())

    return dist

def rms_coordinates(points):

    assert points.ndim == 2

    rms = np.sqrt(
        np.power(points,2).sum(axis=1).mean()
        )

    return rms

def rmsd_coordinates(points_1, points_2):
    """Calculate the RMSD between two sets of coordinates."""

    p1 = np.array(points_1)
    p2 = np.array(points_2)

    assert p1.shape == p2.shape

    diffs = (p1 - p2)

    # NOTE: need to sum over first axis before averaging!

    rmsd = rms_coordinates(diffs)

    return rmsd

def rmsf_coordinates(points):
    """Calculate the RMSF of a set of points around the mean"""

    return rms_coordinates(
        points = (
            points - points.mean(axis=0)
            ),
        )

def pairwise_distances(points_1, points_2=None):
    """Calculate pairwise distances between two sets of points"""

    p1 = (
        np.array(points_1)
        )

    p2 = (
        np.array(points_2)
        if
        points_2 is not None
        else
        np.array(points_1)
        )

    assert p1.ndim == 2
    assert p2.ndim == 2

    assert p1.shape[1] == p2.shape[1], 'not compatible spatial dimensions'

    p1_exp = p1.reshape(
        p1.shape[0], 1, p1.shape[1],
        ).repeat(
        p2.shape[0],
        axis = 1,
        )

    diffs = (p1_exp - p2)

    dists = np.sqrt(
        np.sum(
            np.power(diffs, 2),
            axis=2,
            )
        )

    return dists
