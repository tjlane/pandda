class CommonSetMillerArrayTruncator(object):

    def __init__(self, miller_arrays):

        self.common_set = self.get_miller_arrays_common_set(
            miller_arrays = miller_arrays,
            )

    def __call__(self, miller_array):

        miller_array_truncated = miller_array.common_set(
            self.common_set, 
            assert_is_similar_symmetry = False,
            )

        return miller_array_truncated

    def get_miller_arrays_common_set(self, miller_arrays):

        common_set = miller_arrays[0].set()

        if len(miller_arrays) == 1: 
            return common_set

        for ma in miller_arrays[1:]:

            common_set = common_set.common_set(
                ma, 
                assert_is_similar_symmetry = False,
                )

        return common_set


# class TruncateMillerArraysToResolutionRange:

#     def __init__(self,
#         resolution_high,
#         resolution_low,
#         ):
        
#         self.resolution_high = resolution_high
#         self.resolution_low = resolution_low

#     def __call__(self, miller_array):

        





