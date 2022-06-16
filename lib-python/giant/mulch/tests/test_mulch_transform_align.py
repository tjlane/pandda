import giant.logs as lg
logger = lg.getLogger(__name__)

import copy
import pytest

def test_AlignDatasets(five_baz2b_test_datasets_mcd_labelled):

    # Make a copy as we will be modifying the object
    five_baz2b_test_datasets_mcd_labelled = copy.deepcopy(five_baz2b_test_datasets_mcd_labelled)

    reference_dataset = five_baz2b_test_datasets_mcd_labelled.datasets["BAZ2BA-x434.dimple.pdb"]

    from giant.mulch.transform.align import AlignDatasets
    
    ##

    if True:

        align_datasets = AlignDatasets(
            method = "local",
            alignment_kwargs = dict(
                require_hierarchies_identical = False,
                ),
            processor = None,
            )

        alignments = align_datasets(
            mcd = five_baz2b_test_datasets_mcd_labelled,
            reference_dataset = reference_dataset,
            )

        d_alignment = alignments["BAZ2BA-x430.dimple.pdb"]

        assert d_alignment.n_alignments == 1
        assert d_alignment.n_sites == 115
        assert len(d_alignment.local_alignments) == 1
        assert d_alignment.local_alignments[0].n_sites == 115
        assert d_alignment.local_alignments[0].alignment_mxs[14] is d_alignment.alignment_mxs[14]

        assert list(d_alignment.alignment_mxs[14].r) == pytest.approx(
            [0.999874, 0.002482, 0.015690, -0.002876, 0.999680, 0.025137, -0.015623, -0.025179, 0.999561], abs=1e-6,
            )
        assert list(d_alignment.alignment_mxs[14].t) == pytest.approx(
            [-0.579589, -0.544498, 1.165509], abs=1e-6,
            )
        
        test_site = (24.132558207982353, 18.676095002909786, 24.872317348333294)

        assert test_site == pytest.approx(
            tuple(d_alignment.alignment_sites_ref[17])
            )

        assert test_site == pytest.approx(
            tuple(d_alignment.nat2ref([d_alignment.alignment_sites_nat[17]])[0])
            )

        assert test_site == pytest.approx(
            tuple((d_alignment.alignment_mxs[17].r * d_alignment.alignment_sites_nat[17]) + d_alignment.alignment_mxs[17].t)
            )

    ## 

    if True: 

        align_datasets = AlignDatasets(
            method = "global",
            processor = None,
            )

        alignments = align_datasets(
            mcd = five_baz2b_test_datasets_mcd_labelled,
            reference_dataset = reference_dataset,
            )

        d_alignment = alignments["BAZ2BA-x430.dimple.pdb"]

        assert d_alignment.n_sites == 114 # 

        assert list(d_alignment.alignment_mx.r) == pytest.approx(
            [0.999871, 0.005198, 0.015195, -0.005443, 0.999856, 0.016095, -0.015109, -0.016175, 0.999755], abs=1e-6,
            )

        assert list(d_alignment.alignment_mx.t) == pytest.approx(
            [-0.633184, -0.273177, 0.902664], abs=1e-6,
            )

        test_site = (20.491360228210503, 18.047961342033865, 23.610358031767483)

        assert test_site == pytest.approx(
            tuple(d_alignment.alignment_sites_ref[17])
            )

        assert test_site == pytest.approx(
            tuple(d_alignment.nat2ref([d_alignment.alignment_sites_nat[17]])[0])
            )

        assert test_site == pytest.approx(
            tuple((d_alignment.alignment_mx.r * d_alignment.alignment_sites_nat[17]) + d_alignment.alignment_mx.t)
            )

    ##

    # test errors 
    
    if True:

        align_datasets = AlignDatasets(
            method = "global",
            processor = None,
            alignment_kwargs = {'non_existent_kwarg' : 10},
            )

        alignments = align_datasets(
            mcd = five_baz2b_test_datasets_mcd_labelled,
            reference_dataset = reference_dataset,
            )

        assert sorted(alignments.keys()) == [
            "BAZ2BA-x430.dimple.pdb",
            "BAZ2BA-x431.dimple.pdb",
            "BAZ2BA-x432.dimple.pdb",
            "BAZ2BA-x433.dimple.pdb",
            "BAZ2BA-x434.dimple.pdb",
            ]

        assert list(alignments.values()) == [
            "align_structures_rigid() got an unexpected keyword argument 'non_existent_kwarg'"
            ] * 5
