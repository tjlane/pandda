import os, collections
import numpy, pandas

from pandemic.adp.html import (
    divs,
    HtmlSummary,
    )


class AssessHierarchyGroupsHtmlSummary(HtmlSummary):


    def __init__(self, task):
        self.task = task

    def main_summary(self):

        output = divs.Tab(
            title = 'Analysis of the Hierarchical Partitioning',
            alt_title = 'Hierarchy Analysis',
        )

        tab_set = output.append(divs.TabSet())
        tab_set.append(self.input_vs_output_partitions())
        tab_set.set_active()

        return [output]

    def input_vs_output_partitions(self):

        of = self.task.result.output_files

        output = divs.Tab(alt_title='Filtering of Input Hierarchy')

        txt = """
        > Filtering of hierarchy definitions (input vs optimised)
        Any group that is assigned no disorder in the optimised model is removed from the hierarchy definitions.
        These new eff files are written to the output folder and can be used for further runs and may speed model optimisation.
        """
        output.extend(self.format_summary(txt, classes=["square-corners-top"]))

        #
        # .eff files
        #

        txt = """
        Groups filtered where no atom has a B-factor above <strong>{min_b}&#8491;&#178;</strong>.
        Input Definitions: {in_eff}
        Filtered Definitions: {out_eff}
        """.format(
            min_b = self.task.min_b_factor,
            in_eff = of.get('input_eff_file'),
            out_eff = of.get('output_eff_file'),
            )
        output.extend(self.format_summary(txt))

        #
        # Input Partitions
        #

        # Create output block for this section
        output_block = output.append(divs.Block())
        # Add tab set to the block
        tab_set = output_block.append(divs.TabSet())

        # Add tab for each chain
        input_partitions = of.get('input_partitions_png',{})
        output_partitions = of.get('output_partitions_png',{})
        chain_ids = sorted(set(list(input_partitions.keys()) + list(output_partitions.keys())))

        for c in chain_ids:

            p_in = input_partitions.get(c)
            p_out = output_partitions.get(c)
            tab = divs.Tab(
                alt_title = 'Chain {}'.format(c),
                contents = [
                    divs.Block(width=6, title='Input Hierarchy', image=self.image(p_in)),
                    divs.Block(width=6, title='Filtered Hierarchy', image=self.image(p_out)),
                ],
            )
            tab_set.append(tab)
        tab_set.set_active()

        return output
