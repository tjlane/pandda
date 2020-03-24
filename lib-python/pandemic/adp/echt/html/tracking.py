from libtbx import adopt_init_args
from pandemic.adp.html import HtmlSummary, divs


class EchtTrackingHtmlSummary(HtmlSummary):


    def __init__(self,
        tracking_object,
        ):
        adopt_init_args(self, locals())

    def main_summary(self):

        of = self.tracking_object.output_files

        tab = divs.Tab(title = 'Optimisation Summary')

        ###

        block = tab.append(divs.Block())

        txt = """
        > ECHT model properties during optimisation
        These graphs show how the model characteristics of the hierarchical model change during optimisation.
        At the beginning of optimisation, the number of model parameters and/or the B-factors of each groups are heavilily restrained.
        As these restraints are relaxed over subsequent cycles, the complexity of the model increases, and so the sum of amplitudes increases.
        Towards the end of optimisation, these lines should plateau as the model converges.
        """
        txt_block = divs.Block(
            width = 4,
            contents = self.format_summary(txt),
        )
        block.append(txt_block)

        img_block = divs.Block(
            width = 8,
            image = self.image(of.get('amplitudes_lineplot')),
        )
        block.append(img_block)

        ###

        return [tab]

    def short_summary(self):

        of = self.tracking_object.output_files

        table = self.tracking_object.table.dropna(axis='columns', how='all')

        max_cycle = max(table['cycle'])
        table = table[table['cycle'] == max_cycle]
        table = table.set_index('cycle')

        block = divs.Alert(
            title = 'ECHT statistics at end of optimisation',
            width = 6,
            text = 'Data from {}'.format(of['tracking_csv']),
            table = table.round(1)\
                .to_html(index=False, bold_rows=False, classes=['table table-hover nowrap'])\
                .replace('border="1" ', ''),
        )

        return [block]
