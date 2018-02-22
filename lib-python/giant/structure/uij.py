
import numpy

from scitbx.array_family import flex
from scitbx import linalg

def uij_eigenvalues(uij):
    return linalg.eigensystem_real_symmetric(uij).values()

def check_uij_positive_semi_definite(uij, tol=1e-6):
    """Calculate eigenvalues for each uij and check that all are greater than zero (within tolerance)"""
    # Check tolerance negative
    tol = -1.0 * abs(tol)
    # Extract eigenvalues for each atom
    eigenvalues = numpy.apply_along_axis(uij_eigenvalues, 1, uij)
    # Check all greater than zero
    neg = (eigenvalues < tol)

    return neg.sum()


