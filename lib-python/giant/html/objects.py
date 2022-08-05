import copy
import numpy as np

from . import divs

class ProgressBar(divs.Block):

    __slots__ = (
        'data',
    )

    type = 'progressbar'

    def __init__(
        self,
        data,
        width = 12,
        add_counts = True,
        add_percentages = True,
        **kw_args
        ):

        # MUST BE FIRST
        super(divs.Block, self).__init__(**kw_args)

        self.set(
            width = width,
            )

        self.set_data(
            data = data,
            add_counts = add_counts,
            add_percentages = add_percentages,
            )

    def set_data(self, data, add_counts, add_percentages):

        total_width = float(
            np.sum(
                [d['value'] for d in data]
                )
            )

        self.data = []

        for d in data: 
            
            l = d['label']
            v = d['value']

            p = 100.0 * float(v) / total_width

            if (add_counts and add_percentages):
                l = '{}: {} ({}%)'.format(
                    str(l),
                    str(v),
                    int(p),
                    )
            elif add_counts: 
                l = '{} ({})'.format(
                    str(l), 
                    str(v),
                    )
            elif add_percentages: 
                l = '{} ({}%)'.format(
                    str(l), 
                    int(p),
                    )
            else:
                l = str(l)

            # Pass through all attributes of input dict
            d_copy = copy.deepcopy(d)
            d_copy['size'] = p
            d_copy['text'] = l

            self.data.append(
                d_copy
                )

