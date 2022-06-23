from libtbx import adopt_init_args
from pandemic.adp.html import HtmlSummary, divs


class EchtOptimisationParametersHtmlSummary(HtmlSummary):


    def __init__(self, task):
        adopt_init_args(self, locals())

    def main_summary(self):

        tab = divs.Tab(
            id = 'model_optimisation',
            title = 'Model Optimisation',
        )

        output_files = self.task.result.output_files

        if output_files.get('level amplitudes weights'):

            wgt_block = tab.append(divs.Block())

            txt = """
            > Level amplitude optimisation weights
            Optimisation weights are decayed every cycle, causing an exponential decay in the weight value over cycles.
            This leads to a straight line when plotted on a logarithmic scale.
            When the weights reach a minimum value, they are no decayed any further and so will continue horizontally.
            """
            txt_block = divs.Block(
                width = 4,
                contents = self.format_summary(txt),
            )
            wgt_block.append(txt_block)

            scroll_block = divs.ScrollX(
                width = 8,
            )
            wgt_block.append(scroll_block)

            wgt_dict = output_files.get('level amplitudes weights')
            for variable, image_path in sorted(wgt_dict.items()):
                block = divs.Block(
                    title = 'Weight: {}'.format(variable),
                    image = self.image(image_path),
                    width = 8,
                )
                scroll_block.append(block)

        return [tab]


