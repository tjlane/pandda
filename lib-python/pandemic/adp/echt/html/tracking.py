from libtbx import adopt_init_args
from pandemic.adp.html import HtmlSummary


class EchtTrackingHtmlSummary(HtmlSummary):


    def __init__(self,
        tracking_object,
        ):
        adopt_init_args(self, locals())

    def main_summary(self):

        of = self.tracking_object.output_files

        output = {
            'alt_title' : 'Optimisation Summary',
            'title' : 'Optimisation Summary',
            'fancy_title' : True,
            'contents' : [],
            }

        ###

        block = {'width':12, 'contents' : []}
        output['contents'].append(block)

        txt = """
        > ECHT model properties during optimisation
        """
        txt_block = {
            'width':4,
            'contents' : self.format_summary(txt, classes=['text-justify']),
            }
        block['contents'].append(txt_block)

        img_block = {
            'width' : 8,
            'image' : self.image(of.get('amplitudes_lineplot')),
            }
        block['contents'].append(img_block)

        return [output]

    def short_summary(self):

        of = self.tracking_object.output_files

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
                    'text': 'Data from {}'.format(of['tracking_csv']),
                    'table': table.round(1).to_html(index=False, bold_rows=False, classes=['table table-hover nowrap'])\
                               .replace('border="1" ', ''),
                    },
                ],
            }

        return [panel]
