import copy
from libtbx import adopt_init_args
from libtbx.utils import Sorry, Failure

import numpy
from scitbx.array_family import flex


class UijIsotropicMask(object):


    def __init__(self,
            selection,
            ):

        assert selection.nd() == 1

        n = selection.size()

        n_isotropic = sum(selection)
        n_anisotropic = n - n_isotropic

        all_isotropic = (n == n_isotropic)
        all_anisotropic = (n == n_anisotropic)

        adopt_init_args(self, locals())

    def __call__(self,
            uij_array,
            in_place = False,
            ):

        if (self.all_anisotropic is True): 
            if (in_place is True): 
                return uij_array
            else: 
                return copy.copy(uij_array)

        # Is it a numpy array?
        is_numpy_array = True if hasattr(uij_array, 'shape') else False

        if is_numpy_array is True:
            shape = uij_array.shape[:-1]
            assert uij_array.shape[-1] == 6
            uij_array = flex.sym_mat3_double(uij_array.reshape((numpy.product(shape), 6)))
            uij_array.reshape(flex.grid(shape))
        else:
            shape = uij_array.all()

        if (in_place is False) and (not is_numpy_array):
            uij_array = copy.copy(uij_array)

        # Check compatibility of input array
        if len(shape) == 1:
            assert shape == (self.n,)
            selection = self.selection
        elif len(shape) > 1:
            assert shape[-1] == self.n
            l = numpy.product(shape[:-1])
            selection = (flex.double(l,1).matrix_outer_product(self.selection.as_double()) == 1.0)
            assert selection.all() == (l, self.n)
            selection = selection.as_1d()
        else:
            raise Exception('Not implemented')

        uij_array.reshape(flex.grid((uij_array.size(),)))
        u_sel = uij_array.select(selection)
        sel_n = u_sel.size()
        u_sel = u_sel.as_double().as_numpy_array().reshape((sel_n, 6))
        u_iso = u_sel[:,:3].mean(axis=1).reshape((sel_n,1)).repeat(3, axis=1)
        assert u_iso.shape == (sel_n, 3)
        u_iso = flex.sym_mat3_double(numpy.concatenate([u_iso, numpy.zeros_like(u_iso)], axis=1))
        assert u_iso.all() == (sel_n,)
        uij_array.set_selected(selection, u_iso)
        uij_array.reshape(flex.grid(shape))

        # Convert back to numpy array
        if is_numpy_array:
            uij_array = numpy.array(uij_array).reshape(shape+(6,))

        return uij_array

    def __getitem__(self, mask):
        assert len(mask) == len(self.selection)
        if hasattr(mask, 'shape'):
            mask = flex.bool(mask)
        new_selection = self.selection.select(mask)
        return UijIsotropicMask(selection=new_selection)

    def as_fully_isotropic(self):
        selection = flex.bool(self.selection.size(), True)
        return UijIsotropicMask(selection)

    def as_fully_anisotropic(self):
        selection = flex.bool(self.selection.size(), False)
        return UijIsotropicMask(selection)

    @classmethod
    def from_uij_array(cls,
        array,
        axis = (0,2),
        ):
        selection = flex.bool((array == -1).all(axis=axis))
        return cls(selection)


