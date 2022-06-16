import giant.logs as lg 
logger = lg.getLogger(__name__)

import itertools

import pandas as pd

from giant.processors import (
    basic_processor,
    )

from giant.structure.formatting import (
    short_labeller,
    )

from giant.structure.common import (
    GetInterestingResnames,
    )

from giant.structure.b_factors import (
    CalculateSurroundingsBFactorRatio,
    )

from giant.structure.rmsds import (
    CalculatePairedConformerSetsRMSDs,
    )

from giant.structure.occupancy import (
    CalculateResidueGroupOccupancy,
    )

from .edstats import (
    CalculateEdstatsScores,
    )

from .plots import (
    ValidationRadarPlot,
    ValidationDistributionPlots,
    )


class WriteScores(object):

    def __init__(self,
        make_radar_plot = None,
        make_distributions = None,
        ):

        self.make_radar_plot = (
            make_radar_plot
            if make_radar_plot is not None
            else ValidationRadarPlot()
            )

        self.make_distributions = (
            make_distributions
            if make_distributions is not None
            else ValidationDistributionPlots()
            )

    def __call__(self,
        residue_scores_df,
        dir_path,
        ):

        path_dict = {}

        path_dict.update(
            self.write_residue_scores_csv(
                residue_scores_df = residue_scores_df,
                out_path = (
                    dir_path / "residue_scores.csv"
                    ),
                )
            )

        path_dict.update(
            self.make_residue_scores_radar_plots(
                residue_scores_df = residue_scores_df,
                out_prefix = (
                    dir_path / "residue_scores_"
                    ),
                )
            )

        path_dict.update(
            self.make_residue_scores_distributions(
                residue_scores_df = residue_scores_df,
                out_path = (
                    dir_path / "score_distributions.png"
                    ),
                )
            )

        logger(
            self.format_dict(path_dict)
            )

        return path_dict

    def format_dict(self, path_dict, indent=2):

        indent_str = ' ' * int(indent)

        s_ = "{"

        for k, v in path_dict.items():

            if hasattr(v, 'keys'):

                v_fmt = self.format_dict(v, indent=indent+2)

            else:

                v_fmt = str(v)

            s_ += (
                '\n' + indent_str + "{k} : {v}"
                ).format(
                k = k,
                v = v_fmt,
                )

        s_ += ('\n' + indent_str + "}")

        return s_

    def write_residue_scores_csv(self, residue_scores_df, out_path):

        residue_scores_df.to_csv(
            str(out_path)
            )

        return {'csv' : str(out_path)}

    def make_residue_scores_radar_plots(self, residue_scores_df, out_prefix):

        path_dict = {}

        for index, values in residue_scores_df.iterrows():

            if isinstance(index, tuple):
                label = '_'.join(map(str,index))
            else:
                label = str(index)

            values_dict = values.to_dict()
            values_dict.setdefault('label', label)

            img_path = (
                str(out_prefix) + label + '.png'
                )

            path_dict[index] = str(img_path)

            self.make_radar_plot(
                values_dicts = [values_dict],
                out_path = img_path,
                )

        return {'radar_plots' : path_dict}

    def make_residue_scores_distributions(self, residue_scores_df, out_path):

        scores_df = residue_scores_df.reset_index()
        scores_df = scores_df.dropna(how='all', axis=1)

        if 'Dataset' in scores_df.columns:
            scores_df = scores_df.drop('Dataset', axis=1)

        self.make_distributions(
            scores_df,
            out_path = out_path,
            split_on_column = 'Residue',
            )

        return {'score_distributions' : out_path}


class ScoreModelMultiple(object):

    TABLE_COLUMNS = [
        'Dataset',
        ]

    def __init__(self, 
        score_model = None,
        processor = None,
        ):

        self.score_model = (
            score_model
            if score_model is not None
            else ScoreModel()
            )

        self.processor = (
            processor
            if processor is not None
            else basic_processor
            )

    def __call__(self,
        dataset_list,
        reference_datasets = None,
        ):

        if reference_datasets is None:
            reference_datasets = (
                [None] * len(dataset_list)
                )
        elif len(reference_datasets) == 1:
            reference_datasets = (
                reference_datasets * len(dataset_list)
                )
        else:
            assert len(dataset_list) == len(reference_datasets)

        labels = [
            d.tag for d in dataset_list
        ]

        assert len(labels) == len(set(labels)), (
            'datasets must have unique .tags'
            )

        giant_logger = lg.getLogger('giant')
        logger_level = giant_logger.level
        if logger_level != lg.DEBUG:
            giant_logger.setLevel(max(logger_level,lg.WARNING))

        try:
            dataset_dfs = self.processor([
                self.processor.make_wrapper(
                    func = self.score_model,
                    dataset = d,
                    reference_dataset = r,
                    )
                for d,r in zip(
                    dataset_list,
                    reference_datasets,
                    )
                ])
        finally:
            giant_logger.setLevel(logger_level)

        output_df = self.concatenate_tables(
            tables = dataset_dfs,
            labels = labels,
            )

        logger.subheading('Output scores')

        logger(output_df)

        return output_df

    def concatenate_tables(self, tables, labels):

        df_records = []

        TABLE_COLUMNS = (
            self.TABLE_COLUMNS + self.score_model.TABLE_COLUMNS
            )

        for lab, df in zip(labels, tables):

            df = df.reset_index()

            df['Dataset'] = lab

            df_array = df[TABLE_COLUMNS].to_numpy()

            df_records.extend(df_array)

        df = pd.DataFrame(
            df_records,
            columns = TABLE_COLUMNS,
            )

        df = df.set_index(['Dataset','Residue'])

        df = df.sort_index()

        return df


class ScoreModelSingle(object):

    ROOT_COLUMNS = [
        'RSCC',
        'RSZD',
        'RSZO',
        'RSR',
        'RSZO/OCC',
        'Occupancy',
        'B-factor',
        'B-factor Ratio',
        ]

    TABLE_COLUMNS = [
        'Residue',
        #
        'RSCC',
        'RSZD',
        'RSZO',
        'RSR',
        'RSZO/OCC',
        'Occupancy',
        'B-factor',
        'B-factor Ratio',
        #
        'dRSCC',
        'dRSZD',
        'dRSZO',
        'dRSR',
        'dRSZO/OCC',
        'dOccupancy',
        'dB-factor',
        'dB-factor Ratio',
        #
        'Model RMSD',
    ]

    def __init__(self,
        get_interesting_resnames = None,
        ):

        self.get_interesting_resnames = (
            get_interesting_resnames
            if get_interesting_resnames is not None
            else GetInterestingResnames()
            )

        self.calculate_electron_density_scores = (
            CalculateEdstatsScores(parse_logfile=False)
            )

        self.calculate_surroundings_b_factor_ratio = (
            CalculateSurroundingsBFactorRatio()
            )

        self.calculate_residue_group_occupancy = (
            CalculateResidueGroupOccupancy()
            )

        self.calculate_paired_conformer_sets_rmsds = (
            CalculatePairedConformerSetsRMSDs()
            )

    def __call__(self,
        dataset,
        reference_dataset = None,
        ):

        logger.subheading('Scoring Dataset: {}'.format(dataset.tag))

        dataset_dict = self.score_dataset(
            dataset = dataset,
            )

        if (dataset_dict is None):

            raise ValueError('No scores generated for dataset')

        if (reference_dataset is not None):

            reference_dict = self.score_dataset(
                dataset = reference_dataset,
                )

            self._integrate_reference_dict(
                dataset_dict = dataset_dict,
                reference_dict = reference_dict,
                )

            rmsd_dict = self.calculate_intra_model_rmsds(
                dataset_1 = dataset,
                dataset_2 = reference_dataset,
                )

            self._integrate_rmsd_dict(
                dataset_dict = dataset_dict,
                rmsd_dict = rmsd_dict,
                )

        dataset_dict_list = self._unpack_dict(dataset_dict)

        df = pd.DataFrame(
            data = dataset_dict_list,
            columns = self.TABLE_COLUMNS,
            )

        logger(df)

        return df

    def score_dataset(self, dataset):

        out_dict = dict.fromkeys(self.TABLE_COLUMNS)

        hierarchy = self.sanitise_model(dataset.model.hierarchy)

        residue_groups = self.get_focus_residue_groups(hierarchy)

        if residue_groups is None:
            
            return None

        ed_scores = (
            self.calculate_electron_density_scores(
                dataset = dataset,
                )
            if (
                (dataset.model is not None) and
                (dataset.model.filename is not None) and
                (dataset.data is not None) and
                (dataset.data.filename is not None)
                )
            else None
            )

        logger(str(ed_scores))

        scores_dict = {}

        for rg in residue_groups:

            rg_label = short_labeller.format(rg)

            rg_occ = self.calculate_residue_group_occupancy(
                residue_group = rg,
                )

            rg_dict = {
                'Residue' : rg_label,
                'Occupancy' : rg_occ,
            }

            if ed_scores is not None: 

                rg_scores = ed_scores.get(residue_group=rg)

                rg_dict.update({
                    'RSCC' : rg_scores.rscc,
                    'RSZD' : rg_scores.rszd,
                    'RSZO' : rg_scores.rszo,
                    'RSR'  : rg_scores.rsr,
                    })

                if rg_occ > 0.0:
                    rg_dict['RSZO/OCC'] = (
                        rg_scores.rszo / rg_occ
                        )

            # B-factor ratios to surroundings
            rg_b_scores = self.calculate_surroundings_b_factor_ratio(
                selected_atoms = rg.atoms(),
                hierarchy = hierarchy,
            )

            rg_dict.update({
                'B-factor' : rg_b_scores.selected_average_b,
                'B-factor Ratio' : rg_b_scores.b_factor_ratio,
                })

            scores_dict[rg_label] = (
                rg_dict
                )

        return scores_dict

    def calculate_intra_model_rmsds(self, 
        dataset_1, 
        dataset_2, 
        ):

        h1 = self.sanitise_model(dataset_1.model.hierarchy)
        h2 = self.sanitise_model(dataset_2.model.hierarchy)

        rgs_1 = self.get_focus_residue_groups(h1)
        rgs_2 = self.get_focus_residue_groups(h2)

        rgs_2_dict = {
            short_labeller(rg) : rg
            for rg in rgs_2
        }

        rmsds_dict = {}

        for rg_1 in rgs_1: 

            rg_1_lab = short_labeller(rg_1)

            rg_2 = rgs_2_dict.get(rg_1_lab)

            if rg_2 is None: 
                continue

            confs1, confs2, rmsds = list(zip(
                *self.calculate_paired_conformer_sets_rmsds(
                    conformers_1 = rg_1.conformers(),
                    conformers_2 = rg_2.conformers(),
                )
            ))

            rmsds_dict[rg_1_lab] = (
                max(rmsds)
                )

        return rmsds_dict

    def sanitise_model(self, hierarchy):

        hierarchy = hierarchy.deep_copy()

        hierarchy.atoms().set_chemical_element_simple_if_necessary()
        hierarchy.sort_atoms_in_place()

        return hierarchy

    def get_focus_residue_groups(self, hierarchy):

        resnames = self.get_interesting_resnames(hierarchy)

        if len(resnames) == 0:

            logger('No residue names identified for scoring')

            return None

        resnames_set = set(resnames)

        residue_groups = [
            rg
            for rg in hierarchy.residue_groups()
            if resnames_set.intersection(rg.unique_resnames())
            ]

        return residue_groups

    def _integrate_reference_dict(self, 
        dataset_dict, 
        reference_dict,
        ):

        for r_key, r_dict in reference_dict.items():

            d_dict = dataset_dict.get(r_key)

            if d_dict is None: 
                continue

            for kk in self.ROOT_COLUMNS:

                try:
                    d_dict['d'+kk] = (
                        d_dict[kk] - r_dict[kk]
                        )
                except:
                    pass

    def _integrate_rmsd_dict(self, 
        dataset_dict,
        rmsd_dict,
        ):

        for d_key, d_dict in dataset_dict.items():

            rmsd = rmsd_dict.get(d_key)

            if rmsd is None: 
                continue

            d_dict['Model RMSD'] = (
                rmsd
                )

    def _unpack_dict(self, dataset_dict):

        # Change dict to list of dicts
        dataset_dict_list = []
        for rg_key, rg_dict in sorted(dataset_dict.items()):
            rg_dict['Residue'] = rg_key
            dataset_dict_list.append(rg_dict)

        return dataset_dict_list
