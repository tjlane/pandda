
from libtbx.math_utils import iceil, ifloor

def calculate_minimum_redundancy(unsampled_size, max_sample_size):
    """Calculate the redundancy required to down-sample a vector of size unsampled_size to size max_sample_size"""
    min_redundancy = int(1+(unsampled_size-1)/max_sample_size)
    return min_redundancy

def resample_ordered_list_of_values(vals, redundancy=8):
    """resample a list of values with interpolation"""
    # Number of vals given
    num_inp_vals = len(vals)
    # Number of vals to be returned
    num_samp_vals = int(1+(num_inp_vals-1)/redundancy)
    # Sort in descending order
    ordered_vals = sorted(vals, reverse=True)
    sampled_vals = []

    if num_samp_vals==1:
        return [ordered_vals[0]]
    else:
        sample_dist = (num_inp_vals-1)/(num_samp_vals-1)
#        # Make sure it doesn't overrun the end of the array
#        while sample_dist*(num_samp_vals-1) > num_inp_vals-1:
#            sample_dist = 0.99999999*sample_dist

    # Resample points with interpolation
    for i_point in range(num_samp_vals-1):
        sample_index = sample_dist*i_point
        p1 = ifloor(sample_index)
        v1 = ordered_vals[p1]
        p2 = iceil(sample_index)
        v2 = ordered_vals[p2]
        sample_val = interpolate(x=sample_index, p1=p1, v1=v1, p2=p2, v2=v2)
        sampled_vals.append(sample_val)
    # Add the last point
    sampled_vals.append(ordered_vals[-1])

    assert len(sampled_vals) == num_samp_vals

    return sampled_vals

def interpolate(x,p1,v1,p2,v2):
    """Interpolate the values of p1 and p2 to the point at x"""
    if x==p1: return v1
    if x==p2: return v2
    return v2*(x-p1)/(p2-p1) + v1*(p2-x)/(p2-p1)


