import giant.logs as lg
logger = lg.getLogger(__name__)

import traceback


class GetDummyClass(object):
    """
    A getter for a class that can be passed to override automatically assigned internal classes.
    Typically used to nullify an output class getter.
    """

    def __init__(self, *args, **kwargs):

        pass

    def __call__(self, *args, **kwargs):

        return DummyClass(*args, **kwargs)


class DummyClass(object):
    """
    A class that can be passed to override automatically assigned internal classes.
    Typically used to nullify an output class.
    """

    def __init__(self, *args, **kwargs):

        pass

    def __call__(self, *args, **kwargs):

        pass


class DelayedGetter(object):

    def __init__(self, get_func):

        self.get_func = get_func
        self.func = None

    def __call__(self, *args, **kwargs):

        if (self.func is None):
            self.func = self.get_func()

        return self.func(*args, **kwargs)


class TaskWrapper(object):

    verbose = False

    def __init__(self, func, *args, **kwargs):
        self.func = func
        self.args = args
        self.kwargs = kwargs

    def __call__(self):
        try:
            return self.func(*self.args, **self.kwargs)
        except Exception as e:
            if (self.verbose is True):
                logger(traceback.format_exc())
            raise


class IteratorIterator(object):

    def __init__(self, iterators):
        self.iterators = iterators

    def __iter__(self):
        for i in self.iterators:
            for ii in i: 
                yield ii
        
