import giant.logs as lg
logger = lg.getLogger(__name__)

import os, copy

from libtbx.utils import Sorry, Failure

from giant.mulch.dataset import (
    CrystallographicDataset,
    )


class DefaultReferenceDataset(CrystallographicDataset):

    def __init__(self, model, data):

        super(DefaultReferenceDataset, self).__init__(model=model, data=data)

    def copy(self):
        return copy.deepcopy(self)


class DefaultReferenceSelector(object):

    DatasetClass = DefaultReferenceDataset

    def __init__(self, 
        selection_method = 'resolution',
        pdb_path = None, 
        mtz_path = None, 
        filters = None,
        checks = None,
        process_functions = None,
        ):

        if selection_method not in ["rfree", "resolution"]:
            raise ValueError("Invalid selection_method: {}".format(selection_method))

        self.pdb_path = pdb_path
        self.mtz_path = mtz_path

        self.selection_method = selection_method

        from giant.mulch.filters import DatasetFilterGroup
        self.filter_datasets = DatasetFilterGroup(
            filters = filters,
            )

        from giant.mulch.checks import DatasetCheckerGroup
        self.check_dataset = DatasetCheckerGroup(
            checks = checks,
            )

        from giant.mulch.utils import DatasetProcessor
        self.process_dataset = DatasetProcessor(
            functions = process_functions,
            )

    def __call__(self, mcd=None):

        # Use given reference dataset, or select reference dataset
        if (self.pdb_path and self.mtz_path):
            ref_pdb, ref_mtz = self.pdb_path, self.mtz_path
        else:
            ref_pdb, ref_mtz = self.select_reference_dataset(
                mcd = mcd,
                method = self.selection_method,
                )

        # Load the reference dataset
        reference_dataset = self.load_reference_dataset(
            ref_pdb_path = ref_pdb,
            ref_mtz_path = ref_mtz,
            )

        # return output_paths
        return reference_dataset

    def select_reference_dataset(self, mcd, method='resolution'):
        """Select dataset to act as the reference - scaling, aligning etc"""

        filtered_mcd = self.filter_datasets(mcd)

        if filtered_mcd.n_datasets() == 0: 
            raise ValueError(
                "Can't select a reference dataset - no suitable (non-filtered) datasets remaining"
                )

        if method == 'rfree':
            from giant.mulch.finders import LowestRfreeFinder
            selector = LowestRfreeFinder()
        elif method == 'resolution': 
            from giant.mulch.finders import HighestResolutionFinder
            selector = HighestResolutionFinder()

        i_reference = selector(filtered_mcd)
        reference = mcd.datasets[i_reference]

        return reference.model.filename, reference.data.filename

    def load_reference_dataset(self, ref_pdb_path, ref_mtz_path):
        """Set the reference dataset, to which all other datasets will be aligned and scaled"""

        ref_dataset = self.DatasetClass.from_file(
            model_filename = str(ref_pdb_path),
            data_filename = str(ref_mtz_path),
            ).label(
            num = -1, 
            tag = 'reference',
            )

        # Check selected dataset 
        self.check_dataset(ref_dataset)

        # Apply provided functions
        self.process_dataset(ref_dataset)

        return ref_dataset
