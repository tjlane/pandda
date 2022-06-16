import giant.logs as lg
logger = lg.getLogger(__name__)


class DatasetChecker(object):

    def __init__(self):

        self.name = "DatasetChecker"

    def __call__(self, dataset):
        
        return []


class PossibleStructureFactorChecker(DatasetChecker):

    def __init__(self, possible_structure_factor_pairs):

        super(PossibleStructureFactorChecker, self).__init__()

        self.name = "PossibleStructureFactorChecker"
        self.possible_structure_factor_pairs = possible_structure_factor_pairs

    def __call__(self, dataset):

        mtz_obj = dataset.data.mtz_object()
        dataset_sfs = []
        for sf_pair in self.possible_structure_factor_pairs:
            if mtz_obj.has_column(sf_pair[0]) and mtz_obj.has_column(sf_pair[1]):
                dataset_sfs.append(sf_pair)

        messages = []

        if len(dataset_sfs) > 1:
            messages.append("Dataset contains more than one of the possible structure factors")
        elif len(dataset_sfs) == 0:
            messages.append("Dataset does not contain any of the possible structure factors")

        return messages


class RequiredStructureFactorChecker(DatasetChecker):

    def __init__(self, required_structure_factor_pairs):

        super(RequiredStructureFactorChecker, self).__init__()

        self.name = "RequiredStructureFactorChecker"
        self.required_structure_factor_pairs = required_structure_factor_pairs

    def __call__(self, dataset):

        mtz_obj = dataset.data.mtz_object()
        missing_sfs = []
        for sf_pair in self.required_structure_factor_pairs:
            if not (mtz_obj.has_column(sf_pair[0]) and mtz_obj.has_column(sf_pair[1])):
                missing_sfs.append(sf_pair)

        messages = []
        for miss_pair in missing_sfs:
            messages.append("Dataset is missing structures factors: {}".format(str(miss_pair)))
        return messages


class CrystallographicDatasetModelChecker(DatasetChecker):

    def __init__(self):

        super(CrystallographicDatasetModelChecker, self).__init__()

        self.name = "CrystallographicDatasetModelChecker"

    def __call__(self, dataset):

        messages = []

        if dataset.model.crystal.crystal_symmetry is None:
            messages.append(
                'Could not load crystal symmetry information - check pdb CRYST line is present and valid.'
            )

        if dataset.model.crystal.unit_cell is None:
            messages.append(
                'Could not load unit cell information - check pdb CRYST line is present and valid.'
            )
            
        return messages


class DiffractionDataCompletenessAndValidityChecker(DatasetChecker):

    def __init__(self,
        possible_structure_factor_pairs,
        check_for_invalid_values = True,
        check_completeness_until = 4.0,
        ):

        super(DiffractionDataCompletenessAndValidityChecker, self).__init__()

        self.name = "DiffractionDataCompletenessAndValidityChecker"

        self.possible_structure_factor_pairs = possible_structure_factor_pairs
        self.check_for_invalid_values = check_for_invalid_values
        self.check_completeness_until = check_completeness_until

    def __call__(self, dataset):

        from scitbx.array_family import flex

        mtz_obj = dataset.data.mtz_object()

        dataset_sfs_found = False

        messages = []

        # Find the structure factor from the options
        for sf_pair in self.possible_structure_factor_pairs:

            # Check this is in the dataset
            if mtz_obj.has_column(sf_pair[0]) and mtz_obj.has_column(sf_pair[1]):
                dataset_sfs_found = True
            else: 
                continue

            for col_lab in sf_pair:
                col = mtz_obj.get_column(col_lab)
                ms_d = col.mtz_crystal().miller_set(anomalous_flag=False)
                # Get the boolean selection for valid reflections from the column
                values, valid = col.extract_values_and_selection_valid(-1).as_tuple()
                # List of reflections that have missing values
                zero_sel = flex.bool(valid.size(), True).set_selected(valid, False)

                if (self.check_for_invalid_values is True) and (sum(zero_sel) > 0):

                    # Create list of reflections with zero-values
                    zero_refl = '\n\t'.join([
                        'Reflection: ({:3d}, {:3d}, {:3d}) - resolution: {:6.3f}A - value: {}'.format(
                            i[0], i[1], i[2], d, v
                            )
                        for (i, d), v in
                        list(zip(
                            ms_d.select(zero_sel).d_spacings(),
                            values.select(zero_sel)
                            ))[:10]
                        ])

                    messages.append(
                        'Structure factor column "{}"" has invalid reflections. \n'.format(col_lab) + \
                        '{} reflection(s) are set to N/A or zero (first 10 shown below). \n'.format(sum(zero_sel)) + \
                        'You should populate the structure factors for these reflections with their estimated values\n' + \
                        '(this is normally performed automatically in refinement with phenix or refmac). \n' + \
                        'Analysing maps with missing reflections (especially low resolution reflections!) will degrade the quality of the analysis. \n' + \
                        'However, you can continue by setting checks.all_data_are_valid_values=None. \n' + \
                        'Missing reflections (-1 indicates reflection set to N/A): \n\t{}'.format(zero_refl)
                        )

                if (self.check_completeness_until is not None):

                    # Find selection for the low resolution reflections
                    sel_low = ms_d.resolution_filter_selection(
                        d_min = self.check_completeness_until,
                        d_max = 999.,
                        )
                    # Select these reflections in the dataset miller set
                    ms_d_low = ms_d.select(sel_low)
                    # Extract a complete set of miller indices up to cutoff resolution
                    ms_c_low = ms_d.complete_set(
                        d_min_tolerance = 0.0,
                        d_min = self.check_completeness_until,
                        d_max = 999.)
                    # Calculate the lone set of the two sets
                    lone_set = ms_c_low.lone_set(ms_d_low)

                    # Check that there are the same number of reflections in this set as the other set
                    if lone_set.size() != 0:
                        miss_refl = '\n\t'.join([
                            'Reflection: ({:3d}, {:3d}, {:3d}) - resolution: {:6.3f}A'.format(
                                i[0], i[1], i[2], d
                                )
                             for (i, d)
                             in lone_set.d_spacings()[:10]
                             ])

                        messages.append(
                            'Structure factor column "{}" has missing reflections up to {}A (the column is not complete). \n'.format(col_lab, self.check_completeness_until) + \
                            '{} reflection(s) are missing from the reflection file (first 10 shown below). \n'.format(lone_set.size()) + \
                            'To pass this check you most likely want to expand the reflection set to include unobserved reflections, \n' + \
                            'and then populate the structure factors for these reflections with their estimated values \n' + \
                            '(this is normally performed automatically for 2FOFC/FOFC maps from refinement with phenix or refmac). \n' + \
                            'Analysing maps with missing reflections (especially low resolution reflections!) will degrade the quality of the analysis. \n' + \
                            'However, you can continue by setting all_data_are_valid_values to None. \n' + \
                            '\t{}'.format(miss_refl)
                            )

                    # Calculate overlap between low resolution set and valid set to ensure none missing from low resolution set
                    zero_sel_low = zero_sel.select(sel_low)
                    # Check if any low resolution reflections are invalid
                    if sum(zero_sel_low) > 0:
                        zero_refl = '\n\t'.join([
                            'Reflection: ({:3d}, {:3d}, {:3d}) - resolution: {:6.3f}A - value: {}'.format(
                                i[0], i[1], i[2], d, v
                                )
                            for (i, d), v
                            in list(zip(
                                ms_d_low.select(zero_sel_low).d_spacings(),
                                values.select(sel_low).select(zero_sel_low)
                                ))[:10]
                            ])

                        messages.append(
                            'Structure factor column "{}" has invalid reflections below {}A. \n'.format(col_lab, self.check_completeness_until) + \
                            '{} reflection(s) are set to N/A or zero (first 10 shown below). \n'.format(sum(zero_sel_low)) + \
                            'You should populate the structure factors for these reflections with their estimated values\n' + \
                            '(this is normally performed automatically in refinement with phenix or refmac). \n' + \
                            'Analysing maps with missing reflections (especially low resolution reflections!) will degrade the quality of the analysis. \n' + \
                            'You can turn off this check in the input parameters -- though I really would not -- by setting low_resolution_completeness to None. \n' + \
                            'Missing reflections (-1 indicates reflection set to N/A): \n\t{}'.format(zero_refl)
                            )

        if (dataset_sfs_found is False):
            messages.append(
                (
                    "None of the possible_structure_factor_pairs were found in the dataset: {search}"
                    "\nColumns in dataset: {present}"
                    ).format(
                    search = self.possible_structure_factor_pairs,
                    present = mtz_obj.column_labels(),
                    )
                )

        return messages


################


class DatasetCheckerGroup(object):

    def __init__(self, checks=None):

        if (checks is None): 
            checks = []

        self.name = "DatasetCheckerGroup"
        self.checks = checks

    def __call__(self, dataset):

        messages = []
        for c in self.checks:
            c_msgs = c(dataset=dataset)
            messages.extend(c_msgs)

        return messages


class CrystallographicDatasetChecker(DatasetCheckerGroup):

    def __init__(self):

        super(CrystallographicDatasetChecker, self).__init__()

        self.name = "CrystallographicDatasetChecker"

        self.checks.append(
            CrystallographicDatasetModelChecker()
            )


class DiffractionArrayChecker(DatasetCheckerGroup):

    def __init__(self,
        possible_structure_factor_pairs = None,
        required_structure_factor_pairs = None,
        check_completeness_until = 4.0,
        check_for_invalid_values = True,
        ):

        super(DiffractionArrayChecker, self).__init__()

        self.name = "DiffractionArrayChecker"

        assert (possible_structure_factor_pairs is not None) or (required_structure_factor_pairs is not None)

        if (possible_structure_factor_pairs is not None):
            self.checks.append(
                PossibleStructureFactorChecker(
                    possible_structure_factor_pairs = possible_structure_factor_pairs,
                    )
                )

        if (required_structure_factor_pairs is not None):
            self.checks.append(
                RequiredStructureFactorChecker(
                    required_structure_factor_pairs = required_structure_factor_pairs,
                    )
                )

        if (check_completeness_until is not None) or (check_for_invalid_values is True):

            if (possible_structure_factor_pairs is not None):
                self.checks.append(
                    DiffractionDataCompletenessAndValidityChecker(
                        possible_structure_factor_pairs = possible_structure_factor_pairs,
                        check_for_invalid_values = check_for_invalid_values,
                        check_completeness_until = check_completeness_until,
                        )
                    )

            if (required_structure_factor_pairs is not None):
                for sf_pair in required_structure_factor_pairs:
                    self.checks.append(
                        DiffractionDataCompletenessAndValidityChecker(
                            possible_structure_factor_pairs = [sf_pair],
                            check_for_invalid_values = check_for_invalid_values,
                            check_completeness_until = check_completeness_until,
                            )
                    )


################


class DatasetCheckingFailure(Exception):
    pass


class MultiDatasetChecker(object):

    name = "MultiDatasetChecker"

    def __init__(self, dataset_checker, action='raise'):

        action = str(action).lower()
        assert action in ['raise','warn','none',None]

        self.dataset_checker = dataset_checker
        self.action = action
        self.messages = {}

    def __call__(self, mcd):

        for dtag, dataset in mcd.datasets.items():

            d_messages = self.dataset_checker(dataset)

            if len(d_messages) > 0:
                self.messages.setdefault(dtag,[]).extend(d_messages)

        if len(self.messages) > 0: 
            
            if self.action == 'raise':
                raise DatasetCheckingFailure(
                    self.format_messages(self.messages)
                    )
            elif self.action == 'warn':
                logger.warning(
                    self.format_messages(self.messages)
                    )
            else:
                pass

        return self.messages

    def __str__(self):
        return (
            'Type:{name}\n'
            'Messages:\n{messages}\n'
            ).format(
            name = self.name,
            messages = self.format_messages(self.messages),
            )

    def as_filter(self):

        from giant.mulch.filters import ManualDatasetFilter
        return ManualDatasetFilter(rejections=self.messages)

    def format_messages(self, message_dict):
        strings = [
            "Dataset '{k}': \n\t{s}".format(
                k = k, 
                s = '\n'.join(v).replace('\n','\n\t')
                ) 
            for k,v in message_dict.items()
            ]
        return '\n'.join(strings)

