import os
from libtbx import adopt_init_args
from pandemic.adp.html import HtmlSummary


class TrackingHtmlSummary(HtmlSummary):


    def __init__(self,
        tracking_object,
        ):
        adopt_init_args(self, locals())

    def short_summary(self):

        output = [
            {
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
                },
            ]

        table = self.tracking_object.table_by_dataset.dropna(axis='columns', how='all')
        max_cycle = max(table['cycle'])
        table = table[table['cycle'] == max_cycle]
        table = table.set_index('cycle')

        output.append({
                'type'       : 'panel',
                'title'      : 'Model fit RMSDs by dataset at end of optimisation',
                'width'      : 6,
                'contents'   : [
                    {
                        'text': 'All data can be found in {}'.format(self.tracking_object.tracking_csv2),
                        'table': table.round(1).to_html(index=False, bold_rows=False, classes=['table table-striped table-hover nowrap'])\
                               .replace('border="1" ', ''),
                        },
                    ],
                })

        return output


