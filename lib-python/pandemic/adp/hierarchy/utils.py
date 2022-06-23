from libtbx import adopt_init_args

import numpy
from scitbx.array_family import flex


class StructureFactory(object):


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


class MaskedStructureFactory(object):


    def __init__(self,
        master_h, 
        mask,
        ):
        mask = flex.bool(mask)
        sf = StructureFactory(master_h=master_h.select(mask))
        adopt_init_args(self, locals())

    def custom_copy(self, uij=None, iso=None, blank_copy=False):
        return self.sf.custom_copy(uij=uij, iso=iso, mask=None, blank_copy=blank_copy)


class PartitionBordersFactory(StructureFactory):


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


