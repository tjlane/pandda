import os
from libtbx import adopt_init_args
from pandemic.adp.html import HtmlSummary


class TrackingHtmlSummary(HtmlSummary):


    def __init__(self,
        tracking_object,
        ):
        adopt_init_args(self, locals())

    def main_summary(self):

        of = self.tracking_object.output_files

        output = {
            'alt_title' : 'Optimisation Summary',
            'title' : 'Model Optimisation',
            'fancy_title' : True,
            'contents' : [],
        }

        ###

        block = {
            'title' : 'Model fit during optimisation',
            'width' : 8,
            'image' : self.image(of['rmsds_convergence']),
            }
        output['contents'].append(block)

        ###

        chain_block = {'width':12, 'contents' : []}
        output['contents'].append(chain_block)

        txt = """
        > Model changes during last cycle
        """
        txt_block = {
            'width':4,
            'contents' : self.format_summary(txt, classes=['text-justify']),
            }
        chain_block['contents'].append(txt_block)

        tab_set = {
            'type':'tabs',
            'width' : 8,
            'contents':[],
            }
        chain_block['contents'].append(tab_set)

        for c, p in of.get('model_changes',{}).iteritems():
            tab = {
                'alt_title': 'Chain {}'.format(c),
                'contents' : [{'image':self.image(p)}],
                }
            tab_set['contents'].append(tab)
        if tab_set['contents']:
            tab_set['contents'][0]['active'] = True

        ### TODO ADD OTHER FILES

        return [output]

    def short_summary(self):

        output = []

        of = self.tracking_object.output_files

        output.append({
                'type'       : 'panel',
                'title'      : 'Optimisation Summary',
                'contents'   : [
                    {
                        'width' : 6,
                        'image' : self.image(of['level_convergence']),
                        },
                    {
                        'width' : 6,
                        'image' : self.image(of['snapshots']),
                        },
                    {
                        'width' : 12,
                        'text' : 'Tracking data written to {}'.format(of['tracking_csv1']),
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
                        'text': 'Data from {}'.format(of['tracking_csv2']),
                        'table': table.round(1).to_html(index=False, bold_rows=False, classes=['table table-hover nowrap'])\
                               .replace('border="1" ', ''),
                        },
                    ],
                })

        return output

