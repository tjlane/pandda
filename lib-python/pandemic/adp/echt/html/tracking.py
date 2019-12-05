from libtbx import adopt_init_args
from pandemic.adp.html import HtmlSummary


class EchtTrackingHtmlSummary(HtmlSummary):


    def __init__(self,
        tracking_object,
        ):
        adopt_init_args(self, locals())

    def main_summary(self):

        output = {
            'alt_title' : 'Optimisation Summary',
            'title' : 'Optimisation Summary',
            'fancy_title' : True,
            'contents' : [],
            }

        block = {
            'title' : 'ECHT model traits during optimisation',
            'width' : 6,
            'contents' : [
                {'image' : self.image(self.tracking_object.tracking_png)},
                ],
            }
        output['contents'].append(block)

        return [output]

    def short_summary(self):

        table = self.tracking_object.table.dropna(axis='columns', how='all')

        max_cycle = max(table['cycle'])
        table = table[table['cycle'] == max_cycle]
        table = table.set_index('cycle')

        panel = {
            'type'  : 'alert',
            'title' : 'ECHT statistics at end of optimisation',
            'width' : 6,
            'show'  : True,
            'contents'  : [
                {
                    'text': 'Data from {}'.format(self.tracking_object.tracking_csv),
                    'table': table.round(1).to_html(index=False, bold_rows=False, classes=['table table-hover nowrap'])\
                               .replace('border="1" ', ''),
                    },
                ],
            }

        return [panel]
