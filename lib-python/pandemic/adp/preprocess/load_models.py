import giant.logs as lg
logger = lg.getLogger(__name__)

import os

from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry, Failure

from giant.mulch.dataset import CrystallographicModel, AtomicModel


class ModelLoader(object):


    def __init__(self,
            model_type = 'crystallographic',
            labelling = 'foldername',
            ):
        adopt_init_args(self, locals())

        # Store appropriate model labelling function
        from giant.paths import foldername, filename
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

        # Load input structures
        logger.subheading('Building model list -- {} files'.format(len(pdb_files)))
        logger('Labelling function: {}'.format(self.labelling))
        models = []
        for f in pdb_files:
            # Check file exists!
            if not os.path.exists(f):
                raise IOError('Input file "{}" does not exist!'.format(f))
            # Read pdb file to model object
            try:
                m = self.model_class.from_file(f)
            except Exception as e:
                logger('Error while loading file {}: \n\t{}'.format(f, e))
                raise Sorry('Failed to load pdb file -- have you selected the correct model_type? Current model_type is {}.'.format(self.model_type))
            # Generate label from filename
            l = self.label_func(f)
            # do not allow labels beginning with "."
            if isinstance(l, str) and l.startswith('.'):
                logger('Generated label ({}) begins with "." -- labels beginning with "." are not allowed.'.format(l))
                l = None
            if (not l):
                if (len(pdb_files) == 1) and (self.labelling != 'filename'):
                    logger('No (valid) label created for label function: {}'.format(self.labelling))
                    logger('Trying to label by filename instead')
                    self.labelling = 'filename'
                    from giant.paths import filename
                    self.label_func = filename
                    l = self.label_func(f)
                if not self.is_valid_label(l):
                    raise Sorry('No label/invalid label generated using labelling function "{}"\n\tLabel {}\n\tFile {}\nRename input files or use a different labelling function.'.format(self.labelling, l, f))
            m.label(tag=l)
            models.append(m)

        # Check for duplicate labels
        all_labels = [m.tag for m in models]
        unq_labels = sorted(set(all_labels))
        if len(unq_labels) != len(models):
            # Count the duplicates and report
            counts = [(l, all_labels.count(l)) for l in unq_labels]
            dup_labels = [l for l,c in counts if (c > 1)]
            dup_strs = ['{} (found {} times)'.format(l,c) for l,c in counts if (l in dup_labels)]
            # First report labels and filenames
            label_strs = []
            for m in models:
                if (m.tag not in dup_labels):
                    continue
                label_strs.append('got label "{}" from file "{}"'.format(m.tag, m.filename))
            logger('Duplicate labels found:\n\t{}'.format('\n\t'.join(sorted(label_strs))))
            # Report how many duplicates
            logger('Duplicates:\n\t{}'.format('\n\t'.join(dup_strs)))
            # Raise error
            raise Sorry('Duplicate labels generated for current labelling function. Either rename input files or use different dataset-labelling function.')

        # Sort models for convenience
        models = sorted(models, key=lambda m: m.tag)
        logger('\n{} model(s) loaded:'.format(len(models)))
        for m in models:
            logger('\t"{}" from {}'.format(m.tag, m.filename))

        return models

    def is_valid_label(self, label):

        if (label is None):
            return False

        if not isinstance(label, str):
            return False

        if label.startswith('.'):
            return False

        return True
