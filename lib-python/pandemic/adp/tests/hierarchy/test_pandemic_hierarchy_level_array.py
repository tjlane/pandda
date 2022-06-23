from pytest import raises

import iotbx.pdb

from pandemic.resources.structure_test_snippets import pdb_test_structure_atoms

def test_BuildLevelArrayTask():

    import numpy

    from pandemic.adp.hierarchy.level_array import BuildLevelArrayTask

    hierarchy = iotbx.pdb.hierarchy.input(pdb_string=pdb_test_structure_atoms[0]).hierarchy

    bl = BuildLevelArrayTask(
        overall_selection = None,
        )

    bl.run(
        hierarchy = hierarchy,
        level_group_selection_strings = [
            ['chain A'],
            ['resseq 1911:1920','resseq 1921:1930','resseq 1940:1950','resseq 2000:2010'],
            ['resseq {}'.format(i) for i in range(1900, 1955)],
            ],
        )

    assert bl.result.level_group_array.shape == (3, 428)
    #
    assert (bl.result.level_group_array[0] == 0).all()
    #
    assert set(bl.result.level_group_array[1]) == {-1,0,1,2}
    assert numpy.where(bl.result.level_group_array[1]==-1)[0].tolist() == list(range(0,90)) + list(range(259,335))
    assert numpy.where(bl.result.level_group_array[1]==+0)[0].tolist() == list(range(90,171))
    assert numpy.where(bl.result.level_group_array[1]==+1)[0].tolist() == list(range(171,259))
    assert numpy.where(bl.result.level_group_array[1]==+2)[0].tolist() == list(range(335,428))
    assert numpy.where(bl.result.level_group_array[1]==+3)[0].tolist() == []
    #
    assert set(bl.result.level_group_array[2]) == set(range(0, 51))
    assert numpy.where(bl.result.level_group_array[2]==-1)[0].tolist() == []
    assert numpy.where(bl.result.level_group_array[2]==+0)[0].tolist() == list(range(0,4))
    assert numpy.where(bl.result.level_group_array[2]==+1)[0].tolist() == list(range(4,16))
    assert numpy.where(bl.result.level_group_array[2]==21)[0].tolist() == list(range(171,175))
    assert numpy.where(bl.result.level_group_array[2]==49)[0].tolist() == list(range(412,420))
    assert numpy.where(bl.result.level_group_array[2]==50)[0].tolist() == list(range(420,428))
    assert numpy.where(bl.result.level_group_array[2]==51)[0].tolist() == []

    assert bl.result.overall_atom_mask.shape == (428,)
    assert bl.result.overall_atom_mask.all()

    assert len(bl.result.level_group_selection_strings) == 3
    assert bl.result.level_group_selection_strings[0] == ['chain A']
    assert bl.result.level_group_selection_strings[1] == ['resseq 1911:1920', 'resseq 1921:1930', 'resseq 1940:1950']
    assert bl.result.level_group_selection_strings[2] == ['resseq {}'.format(i) for i in range(1900, 1951)]

    bl = BuildLevelArrayTask(
        overall_selection = 'resseq 1800:1915',
        )

    bl.run(
        hierarchy = hierarchy,
        level_group_selection_strings = [
            ['chain A'],
            ['resseq 1911:1920','resseq 1921:1930','resseq 1940:1950','resseq 2000:2010'],
            ],
        )

    assert bl.result.level_group_array.shape == (2, 133)
    #
    assert (bl.result.level_group_array[0] == 0).all()
    #
    set(bl.result.level_group_array[1]) == {-1,0}
    assert numpy.where(bl.result.level_group_array[1]==-1)[0].tolist() == list(range(0,90))
    assert numpy.where(bl.result.level_group_array[1]==+0)[0].tolist() == list(range(90,133))

    assert bl.result.overall_atom_mask.shape == (428,)
    assert numpy.where(bl.result.overall_atom_mask)[0].tolist() == list(range(0,133))

    assert len(bl.result.level_group_selection_strings) == 2
    assert bl.result.level_group_selection_strings[0] == ['chain A']
    assert bl.result.level_group_selection_strings[1] == ['resseq 1911:1920']

    bl = BuildLevelArrayTask(
        overall_selection = 'resseq 1800:1915',
        )

    with raises(Exception) as e:
        bl.run(
            hierarchy = hierarchy,
            level_group_selection_strings = [
                ['chain A'],
                ['resseq 1920:1930'],
                ],
            )
    assert str(e.value) == "Levels have been created that contain no atoms!\nInvalid Levels: [2]"


def test_BuildLevelArrayTask_duplicate_groups():

    import numpy

    from pandemic.adp.hierarchy.level_array import BuildLevelArrayTask

    hierarchy = iotbx.pdb.hierarchy.input(pdb_string=pdb_test_structure_atoms[0]).hierarchy

    input_selections = [
        ['chain A'],
        ['resseq 1911:1920','resseq 1921'],
        ['resseq 1911:1921'],
        ['resseq 1911:1915','resseq 1916:1920','resseq 1921'],
    ]

    bl = BuildLevelArrayTask(
        overall_selection = None,
        )

    bl.run(
        hierarchy = hierarchy,
        level_group_selection_strings = input_selections,
        )

    assert set(bl.result.level_group_array[0]) == {0}
    assert set(bl.result.level_group_array[1]) == {-1,0,1}
    assert set(bl.result.level_group_array[2]) == {-1,0}
    assert set(bl.result.level_group_array[3]) == {-1,0,1,2}

    assert numpy.where(bl.result.level_group_array[0]==-1)[0].tolist() == []
    assert numpy.where(bl.result.level_group_array[1]==-1)[0].tolist() == list(range(0,90)) + list(range(175,428))
    assert numpy.where(bl.result.level_group_array[2]==-1)[0].tolist() == list(range(0,90)) + list(range(175,428))
    assert numpy.where(bl.result.level_group_array[3]==-1)[0].tolist() == list(range(0,90)) + list(range(175,428))

    assert numpy.where(bl.result.level_group_array[1]==+0)[0].tolist() == list(range(90,171))
    assert numpy.where(bl.result.level_group_array[1]==+1)[0].tolist() == list(range(171,175))

    assert numpy.where(bl.result.level_group_array[2]==+0)[0].tolist() == list(range(90,175))

    assert numpy.where(bl.result.level_group_array[3]==+0)[0].tolist() == list(range(90,133))
    assert numpy.where(bl.result.level_group_array[3]==+1)[0].tolist() == list(range(133,171))
    assert numpy.where(bl.result.level_group_array[3]==+2)[0].tolist() == list(range(171,175))

    # Remove upper groups
    bl = BuildLevelArrayTask(
        overall_selection = None,
        remove_duplicate_groups = 'keep_highest_group',
        )

    bl.run(
        hierarchy = hierarchy,
        level_group_selection_strings = input_selections,
        )

    assert set(bl.result.level_group_array[0]) == {0}
    assert set(bl.result.level_group_array[1]) == {-1,0} #!
    assert set(bl.result.level_group_array[2]) == {-1,0}
    assert set(bl.result.level_group_array[3]) == {-1,0,1,2}

    assert numpy.where(bl.result.level_group_array[0]==-1)[0].tolist() == []
    assert numpy.where(bl.result.level_group_array[1]==-1)[0].tolist() == list(range(0,90)) + list(range(171,428)) #!
    assert numpy.where(bl.result.level_group_array[2]==-1)[0].tolist() == list(range(0,90)) + list(range(175,428))
    assert numpy.where(bl.result.level_group_array[3]==-1)[0].tolist() == list(range(0,90)) + list(range(175,428))

    assert numpy.where(bl.result.level_group_array[1]==+0)[0].tolist() == list(range(90,171))
    #!

    assert numpy.where(bl.result.level_group_array[2]==+0)[0].tolist() == list(range(90,175))

    assert numpy.where(bl.result.level_group_array[3]==+0)[0].tolist() == list(range(90,133))
    assert numpy.where(bl.result.level_group_array[3]==+1)[0].tolist() == list(range(133,171))
    assert numpy.where(bl.result.level_group_array[3]==+2)[0].tolist() == list(range(171,175))

    # Remove lower groups
    bl = BuildLevelArrayTask(
        overall_selection = None,
        remove_duplicate_groups = 'keep_lowest_group',
        )

    bl.run(
        hierarchy = hierarchy,
        level_group_selection_strings = input_selections,
        )

    assert set(bl.result.level_group_array[0]) == {0}
    assert set(bl.result.level_group_array[1]) == {-1,0,1}
    assert set(bl.result.level_group_array[2]) == {-1,0}
    assert set(bl.result.level_group_array[3]) == {-1,0,1} #!

    assert numpy.where(bl.result.level_group_array[0]==-1)[0].tolist() == []
    assert numpy.where(bl.result.level_group_array[1]==-1)[0].tolist() == list(range(0,90)) + list(range(175,428))
    assert numpy.where(bl.result.level_group_array[2]==-1)[0].tolist() == list(range(0,90)) + list(range(175,428))
    assert numpy.where(bl.result.level_group_array[3]==-1)[0].tolist() == list(range(0,90)) + list(range(171,428)) #!

    assert numpy.where(bl.result.level_group_array[1]==+0)[0].tolist() == list(range(90,171))
    assert numpy.where(bl.result.level_group_array[1]==+1)[0].tolist() == list(range(171,175))

    assert numpy.where(bl.result.level_group_array[2]==+0)[0].tolist() == list(range(90,175))

    assert numpy.where(bl.result.level_group_array[3]==+0)[0].tolist() == list(range(90,133))
    assert numpy.where(bl.result.level_group_array[3]==+1)[0].tolist() == list(range(133,171))
    #!

