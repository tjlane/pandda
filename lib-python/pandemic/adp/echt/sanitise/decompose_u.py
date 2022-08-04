from scitbx.matrix import sym, diag
from scitbx.linalg import eigensystem_real_symmetric

def get_negative_u_component(u):
    """
    Get the negative semi-definite component of symmetric matrix u, 
    (optionally setting small negative values to `-minimum_eigenvalue`)
    """

    es = eigensystem_real_symmetric(u)
    vals = es.values()
    vecs = es.vectors().as_scitbx_matrix()

    sel = (vals > -0.0)
    vals.set_selected(sel, 0.0)

    # Convert back to symmetric matrix
    u_out = (vecs.transpose() * diag(diag_elems=vals) * vecs).as_sym_mat3()

    return u_out

def get_positive_u_component(u):
    """
    Get the positive semi-definite component of symmetric matrix u, 
    (optionally setting small positive values to `minimum_eigenvalue`)
    """

    es = eigensystem_real_symmetric(u)
    vals = es.values()
    vecs = es.vectors().as_scitbx_matrix()

    sel = (vals < +0.0)
    vals.set_selected(sel, 0.0)

    # Convert back to symmetric matrix
    u_out = (vecs.transpose() * diag(diag_elems=vals) * vecs).as_sym_mat3()

    return u_out

def decompose_u(u):
    """Split u into contributions with positive eigenvalues and negative eigenvalues"""

    es = eigensystem_real_symmetric(u)
    vals = es.values()
    vecs = es.vectors().as_scitbx_matrix()

    # Create copy of eigenvalues
    pos_eigs = vals
    neg_eigs = vals.deep_copy()

    # zero negative eigenvalues/positive eigenvalues
    pos_eigs.set_selected(pos_eigs<0.0, 0.0)
    neg_eigs.set_selected(neg_eigs>0.0, 0.0)

    pos_comp = (vecs.transpose() * diag(diag_elems=pos_eigs) * vecs).as_sym_mat3()
    neg_comp = (vecs.transpose() * diag(diag_elems=neg_eigs) * vecs).as_sym_mat3()

    return pos_comp, neg_comp
