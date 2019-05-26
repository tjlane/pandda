import copy
import numpy

from scitbx.array_family import flex
from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure

def show_file_dict(self, file_dict, indent=0):
    log = self.log
    s = '  '
    for k, v in file_dict.iteritems():
        if isinstance(v, dict):
            log(s*indent + '> {}'.format(k))
            show_file_dict(self, v, indent+1)
        elif isinstance(v, str):
            log(s*indent + '> {}: {}'.format(k, v))
        else:
            log(s*indent + '> {}'.format(k))
            try:
                for vv in v:
                    log(s*(indent+1)+str(vv))
            except:
                log(s*(indent+1)+str(v))


class StructureFactory:


    def __init__(self,
            master_h,
            ):
        adopt_init_args(self, locals())

    def blank_copy(self):
        h = self.master_h.deep_copy()
        h.atoms().set_uij(flex.sym_mat3_double(h.atoms().size(), [0.0]*6))
        h.atoms().set_b(flex.double(h.atoms().size(), 0.0))
        return h

    def custom_copy(self, uij=None, iso=None, mask=None, blank_copy=False):
        """Copy of the master hierarchy with custom ADPs"""
        if blank_copy is True:
            m_h = self.blank_copy()
        else:
            m_h = self.master_h.deep_copy()
        # Extract atoms
        m_a = m_h.atoms()
        # Extract masked atoms
        if mask is not None:
            assert mask.size() == m_a.size()
            m_a = m_a.select(mask)
        # Apply new uijs
        if uij is not None:
            assert len(uij) == m_a.size()
            m_a.set_uij(flex.sym_mat3_double(uij))
        else: 
            m_a.set_uij(flex.sym_mat3_double(m_a.size(), (-1,-1,-1,-1,-1,-1)))
        # Apply new bs
        if iso is not None:
            assert len(iso) == m_a.size()
            m_a.set_b(flex.double(list(iso)))
        else:
            m_a.set_b(flex.double(m_a.size(), -1))
        return m_h


class PartitionBordersFactory(StructureFactory):


    #def partition_boundaries(self, i_level):
    #    """Find the boundaries between the level partitions"""
    #    l = self.level_array[i_level]
    #    mask = self.atom_mask.tolist()
    #    return self.partition_boundaries_custom(atom_labels=l, mask=mask)

    def partition_boundaries(self, atom_labels, mask=None):
        """Find the boundaries for labelled regions"""

        b = self.blank_copy()

        atom_labels = numpy.array(atom_labels)

        if mask is not None: 
            mask = flex.bool(mask)
            assert mask.size() == b.atoms().size()
        else:
            mask = flex.bool(b.atoms().size(), True)
        assert mask.size() == b.atoms().size()
        assert mask.iselection().size() == len(atom_labels)
        
        for i in numpy.unique(atom_labels):
            # Set the b-factors to the membership of the group
            h = self.blank_copy()
            h.select(mask).atoms().set_b(flex.double((atom_labels==i).tolist()))
            # Where are the boundaries of the group?
            sel_int = numpy.array(h.atoms().extract_b(), dtype=int)
            boundaries = numpy.array(b.atoms()[:-1].extract_b(), dtype=bool) + numpy.array((sel_int[:-1]*(1-sel_int[1:]))+(sel_int[1:]*(1-sel_int[:-1])), dtype=bool)
            b.atoms()[:-1].set_b(flex.double(boundaries.tolist()))
        return b


class UijIsotropicMask:


    def __init__(self,
            selection,
            ):
        assert selection.nd() == 1
        n = selection.size()
        adopt_init_args(self, locals())

    def __call__(self,
            uij_array,
            in_place = False,
            ):

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
        elif len(shape) == 2:
            assert shape[1] == self.n
            l = shape[0]
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

    @classmethod
    def from_uij_array(cls,
        array,
        axis = (0,2),
        ):
        selection = flex.bool((array == -1).all(axis=axis))
        return cls(selection)


class xAtomAndDatasetMask:


    def __init__(self,
            atom_mask,
            dataset_mask,
            ):
        n_atoms = atom_mask.size()
        n_datasets = dataset_mask.size()
        n_selected_atoms = flex.sum(atom_mask.as_int())
        n_selected_datasets = flex.sum(dataset_mask.as_int())
        combined_mask = (dataset_mask.as_double().matrix_outer_product(atom_mask.as_double()) == 1.0)
        adopt_init_args(self, locals())

    def __call__(self,
            array,
            ):
        assert array.all() == (self.n_datasets, self.n_atoms)
        sel = array.as_1d().select(self.combined_mask.as_1d())
        sel.reshape(flex.grid((self.n_selection_datasets, self.n_selected_atoms)))
        return sel
