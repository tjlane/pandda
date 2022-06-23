from pytest import approx, raises
import numpy
from scitbx.array_family import flex

from pandemic.resources.structure_test_snippets import pdb_test_structure_short as pdb_str

def test_structure_factory():

    from iotbx.pdb import hierarchy
    pdb_h = hierarchy.input(pdb_string=pdb_str).hierarchy

    from pandemic.adp.hierarchy.utils import StructureFactory
    sf = StructureFactory(master_h=pdb_h)

    mask = pdb_h.atom_selection_cache().selection('resseq 12 or resseq 14')
    assert list(mask.iselection()) == [5, 6, 7, 8, 9, 10, 11, 12, 13, 22, 23, 24, 25, 26, 27, 28]
    n = mask.iselection().size()
    mask_inv = (mask == False)

    # Input atoms (masked and non-masked)
    in_mask_a = pdb_h.atoms().select(mask)
    in_mask_inv_a = pdb_h.atoms().select(mask_inv)

    # Randomly generate new atoms
    new_vals_b   = flex.double(numpy.random.random(n))
    new_vals_uij = flex.sym_mat3_double([numpy.random.random(6) for i in range(n)])

    zer_vals_b   = flex.double(n,  0.)
    zer_vals_uij = flex.sym_mat3_double(n, ( 0., 0., 0., 0., 0., 0.)).as_double()

    nul_vals_b   = flex.double(n, -1.)
    nul_vals_uij = flex.sym_mat3_double(n, (-1.,-1.,-1.,-1.,-1.,-1.)).as_double()

    # Check blank copies
    blank = sf.blank_copy()
    assert list(blank.atoms().extract_xyz()) == list(pdb_h.atoms().extract_xyz())
    assert blank.atoms().extract_b().all_eq(0.0)
    assert blank.atoms().extract_uij().as_double().all_eq(0.0)

    # B only - not blank
    custom = sf.custom_copy(uij=None, iso=new_vals_b, mask=mask, blank_copy=False)
    custom_mask_a = custom.atoms().select(mask)
    custom_mask_inv_a = custom.atoms().select(mask_inv)
    assert custom_mask_a.extract_xyz().as_double().all_eq(in_mask_a.extract_xyz().as_double())
    assert custom_mask_inv_a.extract_xyz().as_double().all_eq(in_mask_inv_a.extract_xyz().as_double())
    assert custom_mask_a.extract_b().all_eq(new_vals_b)
    assert custom_mask_a.extract_uij().as_double().all_eq(-1.)
    assert custom_mask_inv_a.extract_b().all_eq(in_mask_inv_a.extract_b())
    assert custom_mask_inv_a.extract_uij().as_double().all_eq(in_mask_inv_a.extract_uij().as_double())

    # B only - blank
    custom = sf.custom_copy(uij=None, iso=new_vals_b, mask=mask, blank_copy=True)
    custom_mask_a = custom.atoms().select(mask)
    custom_mask_inv_a = custom.atoms().select(mask_inv)
    assert custom_mask_a.extract_xyz().as_double().all_eq(in_mask_a.extract_xyz().as_double())
    assert custom_mask_inv_a.extract_xyz().as_double().all_eq(in_mask_inv_a.extract_xyz().as_double())
    assert custom_mask_a.extract_b().all_eq(new_vals_b)
    assert custom_mask_a.extract_uij().as_double().all_eq(-1.)
    assert custom_mask_inv_a.extract_b().all_eq(0.0)
    assert custom_mask_inv_a.extract_uij().as_double().all_eq(0.0)

    # U only - not blank
    custom = sf.custom_copy(uij=new_vals_uij, iso=None, mask=mask, blank_copy=False)
    custom_mask_a = custom.atoms().select(mask)
    custom_mask_inv_a = custom.atoms().select(mask_inv)
    assert custom_mask_a.extract_xyz().as_double().all_eq(in_mask_a.extract_xyz().as_double())
    assert custom_mask_inv_a.extract_xyz().as_double().all_eq(in_mask_inv_a.extract_xyz().as_double())
    assert custom_mask_a.extract_b().all_eq(-1.)
    assert custom_mask_a.extract_uij().as_double().all_eq(new_vals_uij.as_double())
    assert custom_mask_inv_a.extract_b().all_eq(in_mask_inv_a.extract_b())
    assert custom_mask_inv_a.extract_uij().as_double().all_eq(in_mask_inv_a.extract_uij().as_double())
    # U only - blank
    custom = sf.custom_copy(uij=new_vals_uij, iso=None, mask=mask, blank_copy=True)
    custom_mask_a = custom.atoms().select(mask)
    custom_mask_inv_a = custom.atoms().select(mask_inv)
    assert custom_mask_a.extract_xyz().as_double().all_eq(in_mask_a.extract_xyz().as_double())
    assert custom_mask_inv_a.extract_xyz().as_double().all_eq(in_mask_inv_a.extract_xyz().as_double())
    assert custom_mask_a.extract_b().all_eq(-1.)
    assert custom_mask_a.extract_uij().as_double().all_eq(new_vals_uij.as_double())
    assert custom_mask_inv_a.extract_b().all_eq(0.0)
    assert custom_mask_inv_a.extract_uij().as_double().all_eq(0.0)

    # U&B only - not blank
    custom = sf.custom_copy(uij=new_vals_uij, iso=new_vals_b, mask=mask, blank_copy=False)
    custom_mask_a = custom.atoms().select(mask)
    custom_mask_inv_a = custom.atoms().select(mask_inv)
    assert custom_mask_a.extract_xyz().as_double().all_eq(in_mask_a.extract_xyz().as_double())
    assert custom_mask_inv_a.extract_xyz().as_double().all_eq(in_mask_inv_a.extract_xyz().as_double())
    assert custom_mask_a.extract_b().all_eq(new_vals_b)
    assert custom_mask_a.extract_uij().as_double().all_eq(new_vals_uij.as_double())
    assert custom_mask_inv_a.extract_b().all_eq(in_mask_inv_a.extract_b())
    assert custom_mask_inv_a.extract_uij().as_double().all_eq(in_mask_inv_a.extract_uij().as_double())
    # U&B only - blank
    custom = sf.custom_copy(uij=new_vals_uij, iso=new_vals_b, mask=mask, blank_copy=True)
    custom_mask_a = custom.atoms().select(mask)
    custom_mask_inv_a = custom.atoms().select(mask_inv)
    assert custom_mask_a.extract_xyz().as_double().all_eq(in_mask_a.extract_xyz().as_double())
    assert custom_mask_inv_a.extract_xyz().as_double().all_eq(in_mask_inv_a.extract_xyz().as_double())
    assert custom_mask_a.extract_b().all_eq(new_vals_b)
    assert custom_mask_a.extract_uij().as_double().all_eq(new_vals_uij.as_double())
    assert custom_mask_inv_a.extract_b().all_eq(0.0)
    assert custom_mask_inv_a.extract_uij().as_double().all_eq(0.0)

    with raises(Exception):
        custom = sf.custom_copy(uij=new_vals_uij, iso=new_vals_b, mask=mask[:-1], blank_copy=True)

    with raises(Exception):
        custom = sf.custom_copy(uij=new_vals_uij, iso=new_vals_b[:-1], mask=mask, blank_copy=True)

    with raises(Exception):
        custom = sf.custom_copy(uij=new_vals_uij[:-1], iso=new_vals_b, mask=mask, blank_copy=True)

def test_partition_borders():

    from iotbx.pdb import hierarchy
    pdb_h = hierarchy.input(pdb_string=pdb_str).hierarchy

    from pandemic.adp.hierarchy.utils import PartitionBordersFactory
    sf = PartitionBordersFactory(master_h=pdb_h)

    atom_labels = [a.resseq for a in pdb_h.atoms_with_labels()]
    bounds = sf.partition_boundaries(atom_labels)
    assert list((bounds.atoms().extract_b() > 0.0).iselection()) == [4, 13, 21, 28]

    atom_labels = [a.resname for a in pdb_h.atoms_with_labels()]
    bounds = sf.partition_boundaries(atom_labels)
    assert list((bounds.atoms().extract_b() > 0.0).iselection()) == [4, 13, 21, 28]

    sel = pdb_h.atom_selection_cache().selection('resseq 12:16')
    assert list(sel.iselection()) == list(range(5,37))
    atom_labels = [a.resseq for a in pdb_h.select(sel).atoms_with_labels()]
    bounds = sf.partition_boundaries(atom_labels, mask=sel)
    assert list((bounds.atoms().extract_b() > 0.0).iselection()) == [4, 13, 21, 28]

    sel = pdb_h.atom_selection_cache().selection('resseq 12 or resseq 14')
    atom_labels = [a.resseq for a in pdb_h.select(sel).atoms_with_labels()]
    bounds = sf.partition_boundaries(atom_labels, mask=sel)
    assert list((bounds.atoms().extract_b() > 0.0).iselection()) == [4, 13, 21, 28]

    sel = pdb_h.atom_selection_cache().selection('resseq 12 or resseq 15')
    atom_labels = [a.resseq for a in pdb_h.select(sel).atoms_with_labels()]
    bounds = sf.partition_boundaries(atom_labels, mask=sel)
    assert list((bounds.atoms().extract_b() > 0.0).iselection()) == [4, 13, 28]

    sel = pdb_h.atom_selection_cache().selection('resseq 12')
    atom_labels = [a.resseq for a in pdb_h.select(sel).atoms_with_labels()]
    bounds = sf.partition_boundaries(atom_labels, mask=sel)
    assert list((bounds.atoms().extract_b() > 0.0).iselection()) == [4, 13]

    sel = pdb_h.atom_selection_cache().selection('resseq 12:13 and name O')
    atom_labels = [a.name for a in pdb_h.select(sel).atoms_with_labels()]
    bounds = sf.partition_boundaries(atom_labels, mask=sel)
    assert list((bounds.atoms().extract_b() > 0.0).iselection()) == [7, 8, 16, 17]
