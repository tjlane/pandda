import giant.logs as lg 
logger = lg.getLogger(__name__)


class AlignDatasetToReference(object):

    def __init__(self, 
        reference_model, 
        method, 
        noerror = True, 
        alignment_kwargs = None,
        ):
        """
        Shortcut object for aligning datasets.
        Constructed so that the class can be initialised and then called within a multiprocessing function.
        """

        if (alignment_kwargs is None):
            alignment_kwargs = {}

        self.reference_model = reference_model
        self.method = method
        self.noerror = noerror
        self.alignment_kwargs = alignment_kwargs

        if method not in ['local', 'global']: 
            raise ValueError(
                'method not defined: {method}'.format(
                    method = self.method,
                    )
                )

    def __call__(self, dataset):

        try: 
            alignment = dataset.model.align_to(
                other_hierarchy = self.reference_model.hierarchy, 
                method = self.method,
                **self.alignment_kwargs
                )
        except Exception as e: 
            if (self.noerror is True): 
                return str(e)
            else: 
                raise

        return alignment


class LabelledAlignDatasetToReference(AlignDatasetToReference): 

    def __call__(self, dataset):

        alignment_or_error = super(
            LabelledAlignDatasetToReference, self,
            ).__call__(
            dataset = dataset,
            )

        return (dataset.tag, alignment_or_error)


class AlignDatasets(object): 

    def __init__(self, method="local", alignment_kwargs=None, processor=None):

        if (processor is None):
            from giant.processors import basic_processor
            processor = basic_processor

        self.method = method 
        self.alignment_kwargs = alignment_kwargs

        self.processor = processor
        self.errors = {}

    def __call__(self, mcd, reference_dataset):

        align_to_reference = LabelledAlignDatasetToReference(
            reference_model = reference_dataset.model,
            method = self.method,
            alignment_kwargs = self.alignment_kwargs,
            )

        arg_list = []
        for dtag, dataset in mcd.datasets.items():

            if not hasattr(dataset, 'tag'):
                dataset.label(tag=dtag)

            wrapper = self.processor.make_wrapper(func=align_to_reference, dataset=dataset)
            arg_list.append(wrapper)

        alignments = self.processor(funcs=arg_list)

        for dtag, align_obj in alignments:

            if isinstance(align_obj, str):
                logger('Dataset {dtag} failed in alignment: \n{msg}'.format(dtag=dtag, msg=align_obj))
                self.errors[dtag] = str(align_obj)
                continue

            dataset = mcd.datasets[dtag]
            dataset.model.alignment = align_obj

        return dict(alignments) 

    def as_filter(self): 

        from giant.mulch.filters import ManualDatasetFilter
        failure_filter = ManualDatasetFilter(
            rejections = self.errors,
            )

        return failure_filter