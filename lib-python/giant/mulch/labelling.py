from giant.paths import foldername, filename


class PathLabeller(object):

    counter_format = '{:05d}'

    def __init__(self, method):

        self.method = method
        
        self.counter = 0

        self.label_func = (
            foldername 
            if (self.method == 'foldername')
            else
            filename
            if (self.method == 'filename')
            else
            self.get_counter
            if (self.method == 'none')
            else
            None
            )

        if self.label_func is None:
            raise NotImplemented()

    def __call__(self, path):

        return self.label_func(
            str(path)
            )

    def get_counter(self, path):

        self.counter += 1

        return self.counter_format.format(
            self.counter
            )