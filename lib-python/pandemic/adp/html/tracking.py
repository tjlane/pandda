import os
from libtbx import adopt_init_args
from pandemic.adp.html import HtmlSummary


class TrackingHtmlSummary(HtmlSummary):


    def __init__(self,
        tracking_object,
        ):
        adopt_init_args(self, locals())

    def main_summary(self):

        output = {
            'alt_title' : 'Optimisation Summary',
            'title' : 'Model Optimisation',
            'fancy_title' : True,
            'contents' : [],
        }

        block = {
            'title' : 'Model fit during optimisation',
            'width' : 6,
            'contents' : [
                {'image' : self.image(self.tracking_object.tracking_png3)},
                ],
            }
        output['contents'].append(block)

        return [output]

    def short_summary(self):

        output = []

        output.append({
                'type'       : 'panel',
                'title'      : 'Optimisation Summary',
                'contents'   : [
                    {
                        'width' : 6,
                        'image' : self.image(
                            self.tracking_object.tracking_png1,
                            ),
                        },
                    {
                        'width' : 6,
                        'image' : self.image(
                            self.tracking_object.tracking_png2,
                            ),
                        },
                    {
                        'width' : 12,
                        'text' : 'Tracking data written to {}'.format(self.tracking_object.tracking_csv1),
                        },
                    ],
                })

        table = self.tracking_object.table_by_dataset.dropna(axis='columns', how='all')
        max_cycle = max(table['cycle'])
        table = table[table['cycle'] == max_cycle]
        table = table.set_index('cycle')

        output.append({
                'type'       : 'alert',
                'title'      : 'Model fit RMSDs by dataset at end of optimisation',
                'width'      : 6,
                'contents'   : [
                    {
                        'text': 'Data from {}'.format(self.tracking_object.tracking_csv2),
                        'table': table.round(1).to_html(index=False, bold_rows=False, classes=['table table-hover nowrap'])\
                               .replace('border="1" ', ''),
                        },
                    ],
                })

        return output

