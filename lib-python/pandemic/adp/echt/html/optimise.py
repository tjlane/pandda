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

            p = divs.Panel(title='Level amplitude optimisation weights')
            tab.append(p)

            scroll = divs.ScrollX()
            p.append(scroll)

            wgt_dict = output_files.get('level amplitudes weights')

            for variable, image_path in sorted(wgt_dict.iteritems()):
                block = divs.Block(
                    title = 'Weight: {}'.format(variable),
                    image = self.image(image_path),
                    width = 5,
                )
                scroll.append(block)

        return [tab]


