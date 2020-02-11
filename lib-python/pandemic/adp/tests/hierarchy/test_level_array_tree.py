
def test_BuildLevelArrayAsTreeTask():

    import numpy

    from pandemic.adp.hierarchy.level_array_tree import BuildLevelArrayAsTreeTask

    level_group_array = -1 * numpy.ones((3, 20))

    level_group_array[0, 00:15] = 0

    level_group_array[1, 00:05] = 0
    level_group_array[1, 10:15] = 1

    level_group_array[2, 00:02] = 0
    level_group_array[2, 02:05] = 1
    level_group_array[2, 05:07] = 2
    level_group_array[2, 07:10] = 3
    level_group_array[2, 10:13] = 4
    level_group_array[2, 13:15] = 5
    level_group_array[2, 15:18] = 6

    bl = BuildLevelArrayAsTreeTask()

    bl.run(level_group_array = level_group_array)

    assert bl.result.tree.links == {
        0: {0: {1: [0, 1], 2: [2, 3]}},
        1: {0: {2: [0, 1]}, 1: {2: [4, 5]}},
        2: {0: {}, 1: {}, 2: {}, 3: {}, 4: {}, 5: {}, 6: {}},
    }

