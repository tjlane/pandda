from builtins import object
import giant.logs as lg
logger = lg.getLogger(__name__)


class DatasetProcessor(object):

    def __init__(self, functions=None):

        if (functions is None):
            functions = []

        self.functions = functions

    def __call__(self, dataset):

        for f in self.functions:
            dataset = f(dataset)

        return dataset