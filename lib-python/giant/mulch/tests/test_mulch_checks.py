import os, copy
import pytest

#####

def test_PossibleStructureFactorChecker(five_baz2b_test_datasets):

    from giant.mulch.checks import PossibleStructureFactorChecker

    success_checker = PossibleStructureFactorChecker(
        possible_structure_factor_pairs = [
            ("NOTF", "NOTSIGF"),
            ("F", "SIGF"),
            ("ALSONOTF", "SIGF"),
            ],
        )

    failure_checker_1 = PossibleStructureFactorChecker(
        possible_structure_factor_pairs = [
            ("NOTF1", "NOTSIGF1"),
            ("NOTF2", "NOTSIGF2"),
            ("F", "NOTSIGF2"),
            ],
        )

    failure_checker_2 = PossibleStructureFactorChecker(
        possible_structure_factor_pairs = [
            ("NOTF1", "NOTSIGF1"),
            ("F", "SIGF"),
            ("2FOFCWT", "PH2FOFCWT"),
            ],
        )

    d0 = five_baz2b_test_datasets[0]

    messages = success_checker(dataset=d0)
    assert len(messages) == 0

    messages = failure_checker_1(dataset=d0)
    assert len(messages) == 1

    messages = failure_checker_2(dataset=d0)
    assert len(messages) == 1

def test_RequiredStructureFactorChecker(five_baz2b_test_datasets):

    from giant.mulch.checks import RequiredStructureFactorChecker

    success_checker = RequiredStructureFactorChecker(
        required_structure_factor_pairs = [
            ("F", "SIGF"),
            ("2FOFCWT", "PH2FOFCWT"),
            ],
        )

    failure_checker_1 = RequiredStructureFactorChecker(
        required_structure_factor_pairs = [
            ("F", "SIGF"),
            ("F", "NOTF"),
            ],
        )

    failure_checker_2 = RequiredStructureFactorChecker(
        required_structure_factor_pairs = [
            ("F", "SIGF"),
            ("F", "NOTF"),
            ("F1", "SIGF1"),
            ],
        )

    d0 = five_baz2b_test_datasets[0]

    messages = success_checker(dataset=d0)
    assert len(messages) == 0

    messages = failure_checker_1(dataset=d0)
    assert len(messages) == 1

    messages = failure_checker_2(dataset=d0)
    assert len(messages) == 2

def test_CrystallographicDatasetModelChecker(five_baz2b_test_datasets):

    from giant.mulch.checks import CrystallographicDatasetModelChecker

    checker = CrystallographicDatasetModelChecker()
    
    d_good = copy.deepcopy(five_baz2b_test_datasets[0])

    assert d_good.model.crystal.crystal_symmetry is not None
    assert d_good.model.crystal.unit_cell is not None

    messages = checker(d_good)
    assert len(messages) == 0

    d_bad = copy.deepcopy(d_good)

    d_bad.model.crystal.crystal_symmetry = None

    messages = checker(d_bad)
    assert len(messages) == 1

    d_bad = copy.deepcopy(d_good)
    d_bad.model.crystal.unit_cell = None

    messages = checker(d_bad)
    assert len(messages) == 1

    d_bad = copy.deepcopy(d_good)
    d_bad.model.crystal.crystal_symmetry = None
    d_bad.model.crystal.unit_cell = None

    messages = checker(d_bad)
    assert len(messages) == 2

def test_DiffractionDataCompletenessAndValidityChecker(five_baz2b_test_datasets):

    from giant.mulch.checks import DiffractionDataCompletenessAndValidityChecker

    success_checker = DiffractionDataCompletenessAndValidityChecker(
        possible_structure_factor_pairs = [
            ("2FOFCWT", "PH2FOFCWT"),
            ("NOTF", "NOTSIGF"),
        ],
        check_for_invalid_values = True,
        check_completeness_until = 4.0,
        )

    failure_checker_1 = DiffractionDataCompletenessAndValidityChecker(
        possible_structure_factor_pairs = [
            ("F", "SIGF"),
            ("NOTF", "NOTSIGF"),
        ],
        check_for_invalid_values = True, # These columns contain invalid values
        check_completeness_until = 4.0, 
        )

    failure_checker_2 = DiffractionDataCompletenessAndValidityChecker(
        possible_structure_factor_pairs = [
            ("2FOFCWT", "PH2FOFCWT"),
        ],
        check_for_invalid_values = True,
        check_completeness_until = 1.5, # should now find missing values
        )

    d = five_baz2b_test_datasets[0]

    messages = success_checker(d)
    assert len(messages) == 0

    messages = failure_checker_1(d)
    assert len(messages) == 4

    messages = failure_checker_2(d)
    assert len(messages) == 2

#####

def test_DatasetCheckerGroup(five_baz2b_test_datasets):

    from giant.mulch.checks import DatasetCheckerGroup, RequiredStructureFactorChecker

    checker = DatasetCheckerGroup(
        checks = [
            RequiredStructureFactorChecker(
                required_structure_factor_pairs = [("NOTF", "SIGNOTF")],
                ),
            RequiredStructureFactorChecker(
                required_structure_factor_pairs = [("NOTG", "SIGNOTG")],
                ),
            ],
        )

    d = five_baz2b_test_datasets[0]

    messages = checker(d)

    assert len(messages) == 2
    assert messages[0] == messages[1].replace("NOTG","NOTF")

def test_CrystallographicDatasetChecker(five_baz2b_test_datasets):

    from giant.mulch.checks import CrystallographicDatasetChecker, CrystallographicDatasetModelChecker

    # At the moment one is just a wrapper around the other
    checker = CrystallographicDatasetChecker()
    checker_1 = CrystallographicDatasetModelChecker()
    
    d = copy.deepcopy(five_baz2b_test_datasets[0])
    d.model.crystal.crystal_symmetry = None
    d.model.crystal.unit_cell = None

    messages = checker(d) 

    messages_1 = checker_1(d)
    assert len(messages_1) == 2

    assert messages == messages_1

def test_DiffractionArrayChecker(five_baz2b_test_datasets):

    from giant.mulch.checks import DiffractionArrayChecker

    checker = DiffractionArrayChecker(
        possible_structure_factor_pairs = [
            ("F","SIGF"),
            ("NOTF","NOTSIGF"),
            ],
        required_structure_factor_pairs = [
            ("2FOFCWT","PH2FOFCWT"),
            ],
            check_completeness_until = 1.5,
            check_for_invalid_values = True,
        )

    d = five_baz2b_test_datasets[0]

    messages = checker(dataset=d)

    assert len(messages) == 8

#####

def test_MultiDatasetChecker(five_baz2b_test_datasets):

    from giant.mulch.collection import MultiCrystalDataset

    datasets = {
        i : copy.deepcopy(d) 
        for i,d in enumerate(five_baz2b_test_datasets)
        }

    datasets[3].model.crystal.crystal_symmetry = None
    datasets[3].model.crystal.unit_cell = None

    mcd = MultiCrystalDataset(datasets=datasets)

    from giant.mulch.checks import CrystallographicDatasetModelChecker
    
    checker = CrystallographicDatasetModelChecker()

    out_msg = (
        "Could not load crystal symmetry information - check pdb CRYST line is present and valid."
        "Could not load unit cell information - check pdb CRYST line is present and valid."
        )

    from giant.mulch.checks import MultiDatasetChecker

    ##

    checker_mcd_raise = MultiDatasetChecker(
        dataset_checker = checker,
        action = "raise",
        )

    with pytest.raises(Exception) as e:
        checker_mcd_raise(mcd)
    assert str(e.typename) == "DatasetCheckingFailure"
    assert str(e.value).replace('\n\t','') == ("Dataset '3': " + out_msg)

    ##

    checker_mcd_warn = MultiDatasetChecker(
        dataset_checker = checker,
        action = "warn",
        )

    messages = checker_mcd_warn(mcd)
    assert list(messages.keys()) == [3]
    assert len(list(messages.values())) == 1
    assert str(''.join(messages[3])).replace('\n\t','') == out_msg

    ##

    checker_mcd_none = MultiDatasetChecker(
        dataset_checker = checker,
        action = "none",
        )

    messages = checker_mcd_none(mcd)
    assert list(messages.keys()) == [3]
    assert len(list(messages.values())) == 1
    assert str(''.join(messages[3])).replace('\n\t','') == out_msg

    ##

    checker_filter = checker_mcd_none.as_filter()

    filter_mcd = checker_filter(mcd)

    assert sorted(filter_mcd.datasets.keys()) == [0,1,2,4]
