import os
from libtbx import adopt_init_args
from pandemic.adp.html import HtmlSummary


class ParameterHtmlSummary(HtmlSummary):


    def __init__(self,
        input_command,
        master_phil,
        running_params,
        parameter_files = None,
        ):
        main_output = {
            'id'        : 'settings',
            'alt_title' : 'Settings',
            'title'     : 'Program Parameters',
            'contents'  : [],
            }
        adopt_init_args(self, locals())
        self.add_top_message()
        self.add_phil_summaries()

    def main_summary(self):
        return [self.main_output]

    def add_top_message(self):

        text = """
        Output written to: {out_dir}.
        """.format(
            out_dir = os.path.abspath(self.running_params.output.out_dir)
            )

        self.main_output['contents'].extend(self.format_summary(text))

    def add_phil_summaries(self):

        phil_str = self.master_phil.format(self.running_params).as_str()
        diff_str = self.master_phil.fetch_diff(source=self.master_phil.format(self.running_params)).as_str()

        panels = []

        panels.append({
            'type'  : 'panel',
            'title' : 'Input command',
            'width' : 12,
            'contents'  : [
                {
                    'width':12,
                    'text':'The following can be copy-pasted into a command shell to rerun the program (assuming the filepaths are correct)',
                    },
                {
                    'width':12,
                    'text':self.wrap_string(string=self.input_command, tag='pre'),
                    },
                ],
            })
        panels.append({
            'type'  : 'panel',
            'title' : 'Parameters different from defaults',
            'width' : 12,
            'contents'  : [
                {
                    'width':12,
                    'text':'The following can be copied to a parameter file (<filename>.eff) and used as input (e.g. <samp> pandemic.adp input.eff</samp>)',
                    },
                {
                    'width':12,
                    'text':self.wrap_string(string=diff_str, tag='pre'),
                    },
                ],
            })
        panels.append({
            'type'  : 'panel',
            'title' : 'All Parameters',
            'width' : 12,
            'show'  : False,
            'contents'  : [
                {
                    'width':12,
                    'text':'The following can be copied to a parameter file (<filename>.eff) and used as input (e.g. <samp> pandemic.adp input.eff</samp>)',
                    },
                {
                    'width':12,
                    'text':self.wrap_string(string=phil_str, tag='pre'),
                    },
                ],
            })

        self.main_output['contents'].extend(panels)

        return


