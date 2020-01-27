from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure
from bamboo.common.logs import Log

from giant.dataset import CrystallographicModel, AtomicModel


class ModelLoader:


    def __init__(self,
            model_type = 'crystallographic',
            labelling = 'foldername',
            verbose = False,
            log = None,
            ):
        adopt_init_args(self, locals())

        # Store appropriate model labelling function
        from bamboo.common.path import foldername, filename
        if labelling == 'foldername':
            self.label_func = foldername
        elif labelling == 'filename':
            self.label_func = filename
        else:
            raise Exception('Invalid labelling function: {}'.format(labelling))

        # Store appropriate model class
        if model_type == "crystallographic":
            self.model_class = CrystallographicModel
        else:
            self.model_class = AtomicModel

    def __call__(self,
            pdb_files,
            ):

        log = self.log

        # Load input structures
        log.subheading('Building model list -- {} files'.format(len(pdb_files)))
        models = []
        for f in pdb_files:
            # Read pdb file to model object
            try: 
                m = self.model_class.from_file(f)
            except Exception as e: 
                self.log('Error while loading file {}: \n\t{}'.format(f, e))
                raise Sorry('Failed to load pdb file -- have you selected the correct model_type? Current model_type is {}.'.format(self.model_type))
            # Generate label from filename
            l = self.label_func(f)
            if not l:
                if len(pdb_files) == 1:
                    log('No label created for label function: {}'.format(self.labelling))
                    log('Trying to label by filename instead')
                    self.labelling = 'filename'
                    from bamboo.common.path import filename
                    self.label_func = filename
                    l = self.label_func(f)
                if not l:
                    raise Sorry('No label generated using labelling function "{}"\n\tLabel {}\n\tFile {}'.format(self.labelling, l, f))
            m.label(tag=l)
            models.append(m)

        # Check for duplicate labels
        all_labels = [m.tag for m in models]
        unq_labels = sorted(set(all_labels))
        if len(unq_labels) != len(models):
            counts = [(l, all_labels.count(l)) for l in unq_labels]
            dups = ['{} (found {} times)'.format(l,c) for l,c in counts if c>1]
            raise Sorry('Duplicate labels generated for models: \n\t{}'.format('\n\t'.join(dups)))

        # Sort models for convenience
        models = sorted(models, key=lambda m: m.tag)
        log('{} models loaded'.format(len(models)))

        return models


