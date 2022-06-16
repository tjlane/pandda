import giant.logs as lg
logger = lg.getLogger(__name__)


class MillerArrayScaler(object):

    def __init__(self, scaling_result):

        self.scaling_result = scaling_result

    def __call__(self, miller_array):

        import copy 
        
        assert miller_array.is_complex_array()

        # Extract intensities
        miller_array_int = miller_array.as_intensity_array()
        # Divide through by amplitudes to get phase component
        miller_array_phs = miller_array * (1.0 / miller_array.as_amplitude_array().data())

        # Need to create a copy to preserve the x-values of the original scaling object
        new_scaling = copy.deepcopy(self.scaling_result)
        new_scaling.new_x_values(
            x_values = miller_array_int.d_star_sq().data(),
        )
        miller_array_scaled_int = miller_array_int.array(
            data = new_scaling.transform(
                miller_array_int.data()
            )
        ).set_observation_type_xray_intensity()

        miller_array_scaled = (
            miller_array_phs * \
            miller_array_scaled_int.as_amplitude_array().data()
            )

        # Check for nan values and set to zero (as phase undefined for such reflections)
        miller_array_scaled.data().set_selected(
            (miller_array_int.data() == 0.0),
            0 + 0j
        )

        assert miller_array_scaled.is_complex_array()

        return miller_array_scaled

    def info(self):
        return self.scaling_result.info


class GetIsotropicMillerArrayScaler(object):

    def __init__(self, reference_miller_array=None):

        self.reference_miller_array = None
        self.scaling_factory = None

        if (reference_miller_array is not None):
            self.set_reference_miller_array(reference_miller_array=reference_miller_array)

    def __call__(self, miller_array):

        assert (self.scaling_factory is not None)

        valid, reason = self.miller_array_is_valid(miller_array)

        if (valid is False):
            raise ValueError(reason)

        scaling_result = self.scaling_factory.calculate_scaling(
            miller_array = miller_array.as_intensity_array(),
        )

        return MillerArrayScaler(
            scaling_result = scaling_result,
            )

    def set_reference_miller_array(self, reference_miller_array):
        self.reference_miller_array = reference_miller_array
        self.scaling_factory = self.get_scaling_factory(reference_miller_array)

    def get_scaling_factory(self, miller_array):
        from giant.xray.scaling import IsotropicBfactorScalingFactory
        factory = IsotropicBfactorScalingFactory(
            reference_miller_array = miller_array.as_intensity_array(),
        )
        return factory

    def miller_array_is_valid(self, miller_array):

        if (miller_array.as_amplitude_array().data().as_numpy_array() == 0).any():
            return False, "contains zero amplitude reflections"

        return True, None


class GetMultiDatasetMillerArrayScalers(object):

    def __init__(self,
        get_miller_array,
        get_scaler,
        ):

        self.name = "MultiDatasetDiffractionScaler"

        self.get_miller_array = get_miller_array
        self.get_scaler = get_scaler

        self.errors = {}

    def __call__(self, mcd, reference_dataset=None):
        """Extract amplitudes and phases for creating map"""

        if (reference_dataset is not None):
            if hasattr(self.get_scaler, "set_reference_miller_array"):
                self.get_scaler.set_reference_miller_array(
                    reference_miller_array = self.get_miller_array(reference_dataset),
                    )

        scalers = {}

        for dtag, dataset in mcd.datasets.items():

            miller_array = self.get_miller_array(dataset)

            try:
                scaler_object = self.get_scaler(miller_array)
            except Exception as e:
                self.errors[dtag] = "error during scaling - {}".format(str(e))
                continue

            scalers[dtag] = scaler_object

        return scalers

    def as_filter(self):

        from giant.mulch.filters import ManualDatasetFilter
        failure_filter = ManualDatasetFilter(rejections=self.errors)

        return failure_filter
