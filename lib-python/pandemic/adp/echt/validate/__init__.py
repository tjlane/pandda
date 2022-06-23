import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy
from scitbx.array_family import flex
from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure


class ValidateEchtModel(object):

    debug = False

    number_of_errors_to_print = 3

    def __init__(self,
            uij_tolerance = 0.01,
            ):

        _uij_fmt_str = ', '.join(6*['{{:.{:d}f}}'.format(int(-1*numpy.log10(uij_tolerance)+1.0))])
        _eig_fmt_str = ', '.join(3*['{{:.{:d}f}}'.format(int(-1*numpy.log10(uij_tolerance)+1.0))])

        adopt_init_args(self, locals())

    def __call__(self,
            model_object,
            reference_hierarchy,
            overall_atom_mask,
            ):
        """Check for negative Uij eigenvalues of output models, etc"""

        overall_atom_mask = flex.bool(overall_atom_mask)
        assert overall_atom_mask.size() == reference_hierarchy.atoms().size()

        reference_atoms = reference_hierarchy.select(overall_atom_mask).atoms()
        assert model_object.n_atoms == reference_atoms.size()

        self.check_output_uijs(
            model_object = model_object,
            reference_atoms = reference_atoms,
            )

    def check_output_uijs(self,
        model_object,
        reference_atoms,
        ):

        from giant.structure.uij import uij_are_positive_semi_definite, sym_mat3_eigenvalues
        from giant.structure.formatting import Labeller

        mo = model_object

        for i_l, u_l in enumerate(mo.uijs()):
            assert u_l.shape == (mo.n_datasets, mo.n_atoms, 6)
            for i_d, u_d in enumerate(u_l):
                uij_valid = uij_are_positive_semi_definite(uij=u_d, tol=self.uij_tolerance)
                if uij_valid.all():
                    continue
                err_msgs = []
                for i_a in numpy.where(numpy.logical_not(uij_valid))[0]:
                    err_msg = 'Atom {:>5d} (Atom Label -> {})'.format(i_a+1, Labeller.format(reference_atoms[i_a]))
                    if mo.all_level_types[i_l] == 'tls':
                        tls_g = None
                        # Find the group that contains this atom - brute force it
                        for i_g, sel_g in enumerate(mo.tls_selections[i_l]):
                            if sel_g[i_a] == True: # Have to do equality comparison because of numpy.bool_ problems
                                tls_g = mo.tls_objects[i_l][i_g]
                                break
                        assert tls_g is not None
                        err_msg += '\n\t\t(in TLS Group {} - {})'.format(i_g+1, tls_g.label)
                    elif mo.all_level_types[i_l] == 'adp':
                        pass
                    else:
                        raise Failure('Unrecognised level type: {}'.format(mo.all_level_types[i_l])) # TODO
                    uij_str = self._uij_fmt_str.format(*tuple(u_d[i_a]))
                    eig_str = self._eig_fmt_str.format(*tuple(sym_mat3_eigenvalues(u_d[i_a])))
                    err_msg += '\n\t\tUij -> ({})\n\t\tEigenvalues -> ({})'.format(uij_str, eig_str)
                    err_msgs.append(err_msg)
                # Report, etc
                if (len(err_msgs) > self.number_of_errors_to_print) and not (self.debug is True):
                    n_hidden = len(err_msgs) - self.number_of_errors_to_print
                    err_msgs = err_msgs[:self.number_of_errors_to_print]
                    err_msgs.append('[...] ({} more similar warning{} not shown)'.format(n_hidden, 's' if n_hidden>1 else ''))
                w = 'Level {} ({}): Uijs for dataset {} are not positive-semi-definite ({} atoms)\n\t{}'.format(i_l+1, mo.all_level_names[i_l], i_d+1, uij_valid.sum(), '\n\t'.join(err_msgs))
                logger.warning(w)

