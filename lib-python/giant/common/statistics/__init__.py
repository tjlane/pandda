import numpy as np

def modified_z_scores(values):
    """Z-scores calculated using deviations from the median rather than the mean"""

    values = np.array(values)

    # Calculate deviations from the median
    medn = np.median(values)
    devs = (values - medn)

    # Calculate median of deviations
    mdev = np.median(np.abs(devs))

    return 0.6745 * devs / mdev

def quartiles(values):

    return np.percentile(values, [25,75])

def iqr(values):

    q1, q3 = quartiles(values)

    return q3-q1

def iqr_outliers(values):

    values = np.array(values)

    q1, q3 = quartiles(array)

    iqr = q3 - q1

    return (
        values < (q1-1.5*iqr)
        ) | (
        values > (q3+1.5*iqr)
        )
