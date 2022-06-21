import numpy as np

def make_bdc_map(
    query_map_data, 
    background_map_data, 
    bdc_value,
    ):
    """Calculate query_map_data - bdc*ref_map_data"""

    assert (query_map_data.shape == background_map_data.shape)

    bdc_map_data = (
        (query_map_data) - 
        (background_map_data * bdc_value)
        )

    return bdc_map_data

def calculate_feature_fraction_correlations(
    query_map_data, 
    background_map_data, 
    feature_mask,
    min_fraction = 0.0, 
    max_fraction = 1.0, 
    increment = 0.01, 
    ):
    """
    Estimate the background corrections of a feature in "query_map", defined by grid point indices "feature_idxs".
    Reference should be an equivalent map differing only by not having the feature present. An optional mask of reference_idxs can be given.
    """

    assert (query_map_data.shape == background_map_data.shape)
    assert (query_map_data.shape == feature_mask.shape)
    assert 1.0 >= max_fraction > min_fraction >= 0.0, (
        'Min ({}) and Max ({}) fraction values are not valid'.format(min_fraction, max_fraction)
        )

    # Background values, global and local
    all_background_values = background_map_data
    feature_background_values = all_background_values[feature_mask]

    # Create list of fractions to subtract
    all_feature_fractions = np.arange(
        min_fraction, 
        max_fraction+increment, 
        increment,
        )

    return_values = []

    # Iterate through different amounts of reference map subtraction to estimate feature correction
    for feature_fraction in all_feature_fractions:

        bdc_value = (1.0 - feature_fraction)

        all_bdc_values = make_bdc_map(
            query_map_data = query_map_data, 
            background_map_data = background_map_data, 
            bdc_value = bdc_value,
            )

        feature_bdc_values = all_bdc_values[feature_mask]

        # Calculate the correlations to the reference maps
        feature_correlation = np.corrcoef(
            feature_background_values,
            feature_bdc_values,   
            )[0,1]

        reference_correlation = np.corrcoef(
            all_background_values,
            all_bdc_values, 
            )[0,1]

        return_values.append(
            (feature_fraction, feature_correlation, reference_correlation)
            )

    return list(zip(*return_values))

def find_maximum_series_discrepancy(labels, series_1, series_2):
    """Calculate the point at which two series are maximally different"""

    assert len(series_1) == len(series_2)
    assert len(series_1) == len(labels)

    diffs = [
        series_1[i] - series_2[i] 
        for i in range(len(series_1))
        ]
        
    return labels[diffs.index(max(diffs))]
