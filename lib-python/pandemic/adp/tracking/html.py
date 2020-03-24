import os
from libtbx import adopt_init_args
from pandemic.adp.html import HtmlSummary, divs


class TrackingHtmlSummary(HtmlSummary):


    def __init__(self,
        tracking_object,
        ):
        adopt_init_args(self, locals())

    def main_summary(self):

        of = self.tracking_object.output_files

        output = divs.Tab(
            title = 'Model Optimisation',
            alt_title = 'Optimisation Summary',
        )

        ###

        block = divs.Block(
            width = 8,
            title = 'Model fit during optimisation',
            image = self.image(of['rmsds_convergence']),
        )
        output.append(block)

        ###

        chain_block = divs.Block(width=12)
        output.append(chain_block)

        txt = """
        > Model changes during last cycle
        """
        txt_block = divs.Block(width=4, contents=self.format_summary(txt, classes=['text-justify']))
        chain_block.append(txt_block)

        tab_set = divs.TabSet(width=8)
        chain_block.append(tab_set)

        for c, p in of.get('model_changes',{}).iteritems():
            tab = divs.Tab(
                alt_title = 'Chain {}'.format(c),
                contents = [divs.Block(image=self.image(p))],
            )
            tab_set.append(tab)
        tab_set.set_active()

        ### TODO ADD OTHER FILES

        return [output]

    def short_summary(self):

        output = []

        of = self.tracking_object.output_files

        panel = divs.Panel(
            title = 'Optimisation Summary',
            contents = [
                divs.Block(width=6, image=self.image(of['level_convergence'])),
                divs.Block(width=6, image=self.image(of['snapshots'])),
                divs.Block(width=12, text='Tracking data written to {}'.format(of['tracking_csv1'])),
            ],
        )
        output.append(panel)

        table = self.tracking_object.table_by_dataset.dropna(axis='columns', how='all')
        max_cycle = max(table['cycle'])
        table = table[table['cycle'] == max_cycle]
        table = table.set_index('cycle')

        table_block = divs.Alert(
            title = 'Model fit RMSDs by dataset at end of optimisation',
            width = 6,
            text = 'Data from {}'.format(of['tracking_csv2']),
            table = table.round(1)\
                .to_html(index=False, bold_rows=False, classes=['table table-hover nowrap'])\
                .replace('border="1" ', ''),
        )
        output.append(table_block)

        return output

