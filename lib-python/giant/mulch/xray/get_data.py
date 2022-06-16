import giant.logs as lg
logger = lg.getLogger(__name__)


class LoadFirstValidMillerArrayFromListOfColumnOptions(object):

    def __init__(self, structure_factor_pairs):

        self.structure_factor_pairs = structure_factor_pairs

    def __call__(self, dataset):

        sf_cols = self.structure_factor_pairs

        # Load diffraction data
        mtz_obj = dataset.data.mtz_object()

        # Iterate through possible structure factor pairs
        dataset_sfs = None
        for sf_pair in sf_cols:
            # Check that the data contains the appropriate column
            if mtz_obj.has_column(sf_pair[0]) and mtz_obj.has_column(sf_pair[1]):
                dataset_sfs = sf_pair
                break

        # Raise error if no columns are identified
        if dataset_sfs is None:
            raise ValueError(
                'No matching structure factors were found in the reflection data for the provided dataset. \n' + \
                'Looking for structure factors: \n\t{}\n'.format('\n\t'.join(map(' and '.join, sf_cols))) + \
                'Structure factors in this dataset: \n\t{}\n'.format('\n\t'.join(mtz_obj.column_labels())) + \
                'You may need to change the input structure_factors selection option.'
            )

        miller_array = dataset.data.get_structure_factors(columns=dataset_sfs)

        return miller_array
