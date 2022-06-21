import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy as np

from .output import (
    GetPanddaStatisticalModelOutputter,
    )


class PanddaStatisticalModel(object):

    name = "PanddaStatisticalModel"

    def __init__(self,
        mu_map_data,
        sigma_adjusted_map_data,
        sigma_uncertainty_calculator,
        sigma_uncertainty_map_selection = None,
        dataset_sigma_uncertainties = None,
        cache_sigma_uncertainties = True,
        ):

        if (dataset_sigma_uncertainties is None):
            dataset_sigma_uncertainties = {}

        self.mu_map_data = mu_map_data
        self.sigma_adjusted_map_data = sigma_adjusted_map_data
        self.sigma_uncertainty_calculator = sigma_uncertainty_calculator
        self.sigma_uncertainty_map_selection = sigma_uncertainty_map_selection
        self.dataset_sigma_uncertainties = dataset_sigma_uncertainties
        self.cache_sigma_uncertainties = cache_sigma_uncertainties

    def __str__(self):

        mu_map_string = (
            'Mean: {mean}\n'
            'Min: {min}\n'
            'Max: {max}\n'
            'Std: {std}\n'
            ).format(
            mean = np.mean(self.mu_map_data),
            min = np.min(self.mu_map_data),
            max = np.max(self.mu_map_data),
            std = np.std(self.mu_map_data),
            )

        sigma_adjusted_string = (
            'Mean: {mean}\n'
            'Min: {min}\n'
            'Max: {max}\n'
            'Std: {std}\n'
            ).format(
            mean = np.mean(self.sigma_adjusted_map_data),
            min = np.min(self.sigma_adjusted_map_data),
            max = np.max(self.sigma_adjusted_map_data),
            std = np.std(self.sigma_adjusted_map_data),
            )

        sigma_uncertainty_string = (
            '\n'.join([
                '{k} : {v}'.format(
                    k=k, v=round(v,3),
                    )
                for k, v in sorted(self.dataset_sigma_uncertainties.items())
                ])
            if self.dataset_sigma_uncertainties is not None
            else "None"
            )

        s_ = (
            'Object: {name}\n'
            '| Mu Map Values:\n'
            '|\t{mu_map_string}\n'
            '| Sigma(Adjusted) Map Values:\n'
            '|\t{sigma_adjusted_string}\n'
            '| Sigma(Uncertainty) Values:\n'
            '|\t{sigma_uncertainty_string}\n'
            '`---->'
            ).format(
            name = self.name,
            mu_map_string = (
                mu_map_string.strip().replace('\n','\n|\t')
                ),
            sigma_adjusted_string = (
                sigma_adjusted_string.strip().replace('\n','\n|\t')
                ),
            sigma_uncertainty_string = (
                sigma_uncertainty_string.strip().replace('\n','\n|\t')
                ),
            )

        return s_.strip()

    def get_sigma_uncertainty(self, map_data=None, dataset_key=None):

        if (self.cache_sigma_uncertainties is True):
            if (dataset_key in self.dataset_sigma_uncertainties):
                return self.dataset_sigma_uncertainties[dataset_key]

        if (self.sigma_uncertainty_map_selection is not None):
            map_data = map_data[self.sigma_uncertainty_map_selection]

        unc = self.sigma_uncertainty_calculator(
            map_data = map_data,
            )

        if (self.cache_sigma_uncertainties is True):
            self.dataset_sigma_uncertainties[dataset_key] = unc

        return unc

    def get_z_map(self,
        dataset_key,
        map_data,
        z_map_type = "adjusted+uncertainty",
        ):

        if z_map_type == "none":
            return map_data - self.mu_map_data
        elif z_map_type == "uncertainty":
            unc = self.get_sigma_uncertainty(
                dataset_key = dataset_key,
                map_data = map_data,
                )
            return (map_data - self.mu_map_data) / unc
        elif z_map_type == "adjusted+uncertainty":
            unc = self.get_sigma_uncertainty(
                dataset_key = dataset_key,
                map_data = map_data,
                )
            norm_map = np.sqrt((self.sigma_adjusted_map_data ** 2.) + (unc ** 2.))
            return (map_data - self.mu_map_data) / norm_map
        else:
            raise ValueError("invalid z_map_type provided: {}".format(z_map_type))

    # def get_event_map(self,
    #     map_data,
    #     bdc,
    #     ):

    #     return map_data - (bdc * self.mu_map_data)


class GetPanddaStatisticalModelFitter(object):

    def __init__(self,
        fit_mu,
        fit_sigma_adjusted,
        get_map_grid,
        load_map_mask_name,
        sigma_uncertainty_map_mask_name,
        processor = None,
        ):

        if (processor is None):
            from giant.processors import basic_processor
            processor = basic_processor

        self.fit_mu = fit_mu
        self.fit_sigma_adjusted = fit_sigma_adjusted

        self.get_map_grid = get_map_grid
        self.load_map_mask_name = load_map_mask_name
        self.sigma_uncertainty_map_mask_name = sigma_uncertainty_map_mask_name

        self.processor = processor

    def __call__(self):

        from pandda.analyse.statistical_model.fit import (
            PanddaStatisticalModelFitter,
            )

        map_grid = self.get_map_grid()
        map_mask = map_grid.masks[self.load_map_mask_name]
        unc_mask = map_grid.masks[self.sigma_uncertainty_map_mask_name]
        unc_mask_reindex, unc_mask_selection = unc_mask.index_on_other(map_mask)

        return PanddaStatisticalModelFitter(
            processor = self.processor,
            fit_mu = self.fit_mu,
            fit_sigma_adjusted = self.fit_sigma_adjusted,
            sigma_uncertainty_map_selection = unc_mask_selection,
            )
