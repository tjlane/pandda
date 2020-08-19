

def wrapper_call(func):
    return func()


class ProcessWrapper(object):

    def __init__(self, func, *args, **kw_args):

        self.func = func
        self.args = args
        self.kw_args = kw_args

    def __call__(self):
        return self.func(*self.args, **self.kw_args)


class BaseProcessor(object):

    def make_wrapper(self, func, *args, **kw_args):
        return ProcessWrapper(func=func, *args, **kw_args)

    def as_dict_processor(self):
        return ProcessorDict(processor=self)


class Processor(BaseProcessor): 

    def __init__(self):
        pass

    def __call__(self, funcs):

        results = []

        for func in funcs:
            result = func()
            results.append(result)

        return results


class ProcessorJoblib(BaseProcessor):

    def __init__(self, 
        n_cpus, 
        verbosity=8, 
        ignore_runtime_warnings=True,
        backend = "multiprocessing",
        ):

        self.n_cpus = n_cpus
        self.verbosity = verbosity
        self.ignore_runtime_warnings = ignore_runtime_warnings
        self.backend = backend

        # self.parallel = None

    # def __enter__(self):

    #     import joblib

    #     self.parallel = joblib.Parallel(
    #         n_jobs = self.n_cpus,
    #         verbose = self.verbosity,
    #         ).__enter__()
    #     return self

    # def __exit__(self, exc_type, exc_val, exc_tb):
    #     self.parallel.__exit__(exc_type, exc_val, exc_tb)

    def __call__(self, funcs):

        import joblib, warnings

        with warnings.catch_warnings():
        
            if self.ignore_runtime_warnings: 
                #warnings.simplefilter('ignore', category=RuntimeWarning)
                warnings.simplefilter("ignore")

            results = joblib.Parallel(
                n_jobs = self.n_cpus,
                verbose = self.verbosity,
                backend = self.backend,
                )(
                joblib.delayed(f)() for f in funcs
                )
        
        return results


class ProcessorDict(BaseProcessor):

    def __init__(self, processor=None):
        """
        Wrapper around Processor that allows the supplying of jobs as a dictionary of jobs. 
        Results for each key, func pair will be returned as key, result dictionary.
        """

        if (processor is None):
            processor = Processor()

        self.processor = processor

    def __call__(self, funcs):

        # unpack
        keys = funcs.keys()
        values = [funcs[k] for k in keys]

        processed = self.processor(values)

        results = {
            k: processed[i]
            for i, k 
            in enumerate(keys)
            }

        return results


# initialise some classes for convenience

basic_processor = Processor()


# class ProcessorDask(BaseProcessor):

#     pass

    # with worker_client(timeout=120, separate_thread=False) as client:
    #                 dataset_maps_futures = client.map(wrapper_run, arg_list)
    #                 dataset_maps = client.gather(dataset_maps_futures)    

