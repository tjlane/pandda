import os
from libtbx import adopt_init_args

from pandemic.adp.html import (
    divs,
    HtmlSummary,
    )


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

        chain_block = output.append(divs.Block(width=12))

        txt = """
        > Model changes during last cycle
        The change in the B-factors of each level over the last cycle of optimisation.
        The B-factor changes are averaged over each residue.
        """
        txt_block = divs.Block(width=4, contents=self.format_summary(txt))
        chain_block.append(txt_block)

        tab_set = chain_block.append(divs.TabSet(width=8))
        for c, p in of.get('model_changes',{}).items():
            tab = divs.Tab(
                alt_title = 'Chain {}'.format(c),
                contents = [divs.Block(image=self.image(p))],
            )
            tab_set.append(tab)
        tab_set.set_active()

        ###

        rmsd_block = output.append(divs.Block(width=12))

        txt = """
        > Fit to input ADPs during optimisation
        The root-mean-squared difference between the input ADPs and the hierarchical model ADPs, on the scale of the B-factor.
        At the beginning of optimisation the number of model parameters is limited, so the fit will be poor (high rmsd).
        As the restriction on the number of model parameters is relaxed, the model fit should improve (lower rmsd).
        For a well-fitting model, a final RMSD of ~1 B-factor or less would be expected.
        """
        txt_block = divs.Block(width=4, contents=self.format_summary(txt))
        rmsd_block.append(txt_block)

        img_block = divs.Block(
            width = 8,
            image = self.image(of['rmsds_convergence']),
        )
        rmsd_block.append(img_block)

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

