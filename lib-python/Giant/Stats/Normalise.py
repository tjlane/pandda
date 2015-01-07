import numpy

def normalise_array_to_z_scores(input_array, element_means, element_stds, binary_mask=None):
    """Takes the values in input_array (numpy_array) and scales them by the elements in element_means, element_stds"""

    assert input_array.ndim == 1

    # Number of data points
    num_data = len(input_array)

    if not binary_mask: binary_mask = [1]*num_data
    assert len(binary_mask) == num_data

    # Create empty array to hold output
    z_array = numpy.empty(num_data)

    # Populate the entries in z_array
    [z_array.put(i, (x - element_means[i])/element_stds[i]) if binary_mask[i]==1 else z_array.put(i, 0) for i,x in enumerate(input_array)]

    assert i+1 == num_data

    return z_array
