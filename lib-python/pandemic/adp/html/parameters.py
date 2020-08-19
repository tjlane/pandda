import os
from libtbx import adopt_init_args

from . import (
    divs,
    HtmlSummary,
    )


class ParameterHtmlSummary(HtmlSummary):


    def __init__(self,
        input_command,
        master_phil,
        running_params,
        ):
        main_output = divs.Tab(
            id = 'settings',
            title = 'Program Parameters',
            alt_title = 'Settings',
        )
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

        self.main_output.extend(self.format_summary(text))

    def add_phil_summaries(self):

        phil_str = self.master_phil.format(self.running_params).as_str()
        diff_str = self.master_phil.fetch_diff(source=self.master_phil.format(self.running_params)).as_str()

        output = []

        output.append(
            divs.Panel(
                title = 'Input command',
                contents = [
                    divs.Block(
                        text = 'The following can be copy-pasted into a command shell to rerun the program (assuming the filepaths are correct)',
                    ),
                    divs.Block(
                        text = self.wrap_string(string=self.input_command, tag='pre'),
                    ),
                ],
            )
        )

        output.append(
            divs.Panel(
                title = 'Parameters different from defaults',
                contents = [
                    divs.Block(
                        text = 'The following can be copied to a parameter file (<filename>.eff) and used as input (e.g. <samp> pandemic.adp input.eff</samp>)',
                    ),
                    divs.Block(
                        text = self.wrap_string(string=diff_str, tag='pre'),
                    ),
                ],
            )
        )

        output.append(
            divs.Panel(
                title = 'All Parameters',
                contents = [
                    divs.Block(
                        text = 'The following can be copied to a parameter file (<filename>.eff) and used as input (e.g. <samp> pandemic.adp input.eff</samp>)',
                    ),
                    divs.Block(
                        text = self.wrap_string(string=phil_str, tag='pre'),
                    ),
                ],
            )
        )

        self.main_output.extend(output)

        return


