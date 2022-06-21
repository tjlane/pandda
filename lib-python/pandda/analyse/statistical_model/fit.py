import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy as np

from scipy import stats
from scipy.stats import kde

from .ospina import estimate_true_underlying_sd


class CalculateMu(object):

    def __init__(self):
        pass

    def __call__(self, map_array):

        # Only one dataset?
        if len(map_array) == 1:
            return map_array[0]

        return map_array.mean(axis=0)


class SigmaUncertaintyCalculator(object):

    def __init__(self, n_values=None, q_cut=1.5):

        self.n_values = n_values
        self.q_cut = q_cut

        if n_values is not None:
            self.set_n(n_values)

    def set_n(self, n_values):

        # Extract the theoretical quantiles that we would expect if these values were from a normal distribution
        self.theoretical_values = (
            stats.norm.ppf(np.linspace(0.,1.,n_values+2), 0., 1.)[1:-1] # chop ends which are inf
            )

        # Select the points in the middle of the distribution
        self.middle_selection = (
            np.abs(self.theoretical_values) < self.q_cut
            )

        self.middle_theoretical_values = (
            self.theoretical_values[self.middle_selection]
            )

    def set_reference_data(self, map_data):
        self.reference_map_data = map_data
        self.set_n(len(map_data))

    def __call__(self, map_data):

        if (self.n_values is None) or (len(map_data) != self.n_values):
            self.set_n(len(map_data))

        # Calculate the difference from the reference values
        actual_values = np.array(map_data - self.reference_map_data)
        actual_values.sort()

        middle_actual_values = actual_values[self.middle_selection]

        # Calculate the slope of the centre of the graph
        map_unc, map_off = np.polyfit(
            x = self.middle_theoretical_values,
            y = middle_actual_values,
            deg = 1,
            )

        return map_unc


class CalculateSigmaUncertainties(object):

    def __init__(self, 
        map_selection = None,
        processor = None, 
        ):

        if (processor is None):
            from giant.processors import basic_processor
            processor = basic_processor

        self.calculator = SigmaUncertaintyCalculator()

        self.map_selection = map_selection

        self.processor = processor

    def __call__(self,
        map_array,
        mu_map_data,
        ):

        if (self.map_selection is not None):
            mu_map_data = mu_map_data[self.map_selection]

        self.calculator.set_reference_data(
            map_data = mu_map_data,
            )

        jobs = []

        for i, map_data in enumerate(map_array):

            if (self.map_selection is not None):
                map_data = map_data[self.map_selection]

            func = self.processor.make_wrapper(
                func = self.calculator,
                map_data = map_data,
                )

            jobs.append(func)

        sigma_uncertainties = self.processor(jobs)

        return sigma_uncertainties


class SigmaAdjustedCalculator(object):

    def __init__(self, minimal_guess=1e-6):
        self.minimal_guess = minimal_guess

    def __call__(self, map_array, mu_map_data, sigma_uncertainties):

        import warnings

        sadj_values = np.empty_like(mu_map_data)

        initial_guess = self.minimal_guess

        with warnings.catch_warnings():

            warnings.simplefilter("ignore")

            for i, map_values in enumerate(map_array.T):

                sadj = estimate_true_underlying_sd(
                    obs_vals = map_values,
                    obs_error = sigma_uncertainties,
                    est_mu = mu_map_data[i],
                    est_sigma = initial_guess,
                    )

                sadj_values[i] = sadj

                # Update guess as shortcut
                initial_guess = max(
                    sadj, self.minimal_guess
                )

        return sadj_values


class CalculateSigmaAdjusted(object):

    def __init__(self, processor, chunk_size=1000):

        self.processor = processor
        self.calculator = SigmaAdjustedCalculator()
        self.chunk_size = chunk_size

    def __call__(self,
        map_array,
        mu_map_data,
        sigma_uncertainties,
        ):

        if len(map_array) == 1: 
            sadj_map_data = np.zeros_like(mu_map_data)
            return sadj_map_data

        job_iterator = self.chunk_iterator(
            map_array = map_array,
            mu_map_data = mu_map_data,
            sigma_uncertainties = sigma_uncertainties,
            )

        result_chunks = self.processor(job_iterator)

        sadj_map_data = np.empty_like(mu_map_data)

        for i, map_chunk in enumerate(result_chunks):

            sadj_map_data[i*self.chunk_size:(i+1)*self.chunk_size] = map_chunk

        return sadj_map_data

    def chunk_iterator(self,
        map_array,
        mu_map_data,
        sigma_uncertainties,
        ):

        for i in range(0, len(mu_map_data), self.chunk_size):

            func = self.processor.make_wrapper(
                func = self.calculator,
                map_array = map_array[:,i:i+self.chunk_size],
                mu_map_data = mu_map_data[i:i+self.chunk_size],
                sigma_uncertainties = sigma_uncertainties,
                )

            yield func


class PanddaStatisticalModelFitter(object):

    def __init__(self,
        fit_mu = True,
        fit_sigma_adjusted = True,
        sigma_uncertainty_map_selection = None,
        processor = None,
        ):

        if (processor is None):
            from giant.processors import basic_processor
            processor = basic_processor

        self.processor = processor

        if (fit_mu is True):
            self.calculate_mu = CalculateMu()
        else:
            self.calculate_mu = None

        self.calculate_sigma_uncertainties = CalculateSigmaUncertainties(
            map_selection = sigma_uncertainty_map_selection,
            processor = processor,
            )

        if (fit_sigma_adjusted is True):
            self.calculate_sigma_adjusted = CalculateSigmaAdjusted(
                processor = processor,
                )
        else:
            self.calculate_sigma_adjusted = None

    def __call__(self, map_manager):

        from . import PanddaStatisticalModel

        map_array = map_manager.map_data
        dataset_keys = map_manager.dataset_keys

        ###

        if (self.calculate_mu is not None):
            logger('\n> Calculating Mu Values\n')
            mu_map_data = self.calculate_mu(map_array)
        else:
            logger('\n> Using Zero Mu Values\n')
            mu_map_data = (
                np.zeros(map_array.shape[1], dtype=float)
                )

        ###

        logger('\n> Calculating Sigma-Uncertainty Values\n')
        sigma_uncertainties = self.calculate_sigma_uncertainties(
            map_array = map_array,
            mu_map_data = mu_map_data,
            )

        ###

        if (self.calculate_sigma_adjusted is not None):
            logger('\n> Calculating Sigma-Adjusted Values\n')
            sigma_adjusted_map_data = self.calculate_sigma_adjusted(
                map_array = map_array,
                mu_map_data = mu_map_data,
                sigma_uncertainties = sigma_uncertainties,
                )
        else:
            logger('\n> Using Zero Sigma-Adjusted Values\n')
            sigma_adjusted_map_data = (
                np.zeros(map_array.shape[1], dtype=float)
                )

        ###

        dataset_sigma_uncertainties = dict(
            zip(dataset_keys, sigma_uncertainties)
            )

        ###

        fitted_model = PanddaStatisticalModel(
            mu_map_data = mu_map_data,
            sigma_adjusted_map_data = sigma_adjusted_map_data,
            sigma_uncertainty_calculator = self.calculate_sigma_uncertainties.calculator,
            sigma_uncertainty_map_selection = self.calculate_sigma_uncertainties.map_selection,
            dataset_sigma_uncertainties = dataset_sigma_uncertainties,
            cache_sigma_uncertainties = True,
            )

        logger.subheading('Fitted Statistical Model')
        logger(str(fitted_model))

        return fitted_model
