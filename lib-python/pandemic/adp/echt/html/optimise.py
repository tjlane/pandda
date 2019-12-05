from libtbx import adopt_init_args
from pandemic.adp.html import HtmlSummary


class EchtOptimisationParametersHtmlSummary(HtmlSummary):


    def __init__(self, task):
        adopt_init_args(self, locals())

    def main_summary(self):

        tab = {'id'        : 'model_optimisation',
               'alt_title' : 'Model Optimisation',
               'title' : 'Model Optimisation',
               'fancy_title' : True,
               'contents': [],
              }

        output_files = self.task.result.output_files

        output = []

        if output_files.get('level amplitudes weights'):

            p = {
                'type'      : 'panel',
                'title'     : 'Level amplitude optimisation weights',
                'width'     : 12,
                'contents'  : [],
                }
            output.append(p)

            scroll = {
                'type' : 'scroll',
                'contents' : [],
            }
            p['contents'].append(scroll)

            wgt_dict = output_files.get('level amplitudes weights')

            for variable, image_path in sorted(wgt_dict.iteritems()):
                block = {
                    'title': 'Weight: {}'.format(variable),
                    'image': self.image(image_path),
                    'width': 5,
                }
                scroll['contents'].append(block)

        tab['contents'].extend(output)

        return [tab]


