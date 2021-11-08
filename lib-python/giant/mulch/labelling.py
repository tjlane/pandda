import pathlib as pl

from giant.paths import foldername, filename

def basename(path):

    p = pl.Path(path)

    while p.suffix:
        p = p.with_name(p.stem)

    return p.stem


class PathLabeller(object):

    counter_format = '{:05d}'

    def __init__(self, method='filename'):

        self.method = method
        
        self.counter = 0

        self.label_func = (
            foldername 
            if (self.method == 'foldername')
            else
            filename
            if (self.method == 'filename')
            else
            basename
            if (self.method == 'basename')
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