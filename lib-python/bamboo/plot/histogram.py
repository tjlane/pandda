
from bamboo.plot import finite

from matplotlib import pyplot

def simple_histogram(filename, data, title, x_lab='x', n_bins=30):
    """Quick histogram function"""

    fig = pyplot.figure()
    pyplot.title(title)
    pyplot.hist(x=finite(data), bins=n_bins)
    pyplot.xlabel(x_lab)
    pyplot.ylabel('Count')
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close(fig)

    return

