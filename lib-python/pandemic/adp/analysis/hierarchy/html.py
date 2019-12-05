import os, collections
import numpy, pandas
from libtbx import adopt_init_args
from pandemic.adp.html import HtmlSummary


class AssessHierarchyGroupsHtmlSummary(HtmlSummary):


    def __init__(self, task):
        adopt_init_args(self, locals())

    def main_summary(self):

        output = {
            'alt_title' : 'Hierarchy Analysis',
            'title' : 'Analysis of the Hierarchical Partitioning',
            'fancy_title' : True,
            'contents' : [],
        }

        tab_set = {'type':'tabs', 'contents':[]}
        output['contents'].append(tab_set)

        tab_set['contents'].append(self.input_vs_output_partitions())

        if tab_set['contents']:
            tab_set['contents'][0]['active'] = True
        
        return [output]

    def input_vs_output_partitions(self):

        of = self.task.result.output_files

        output = {
            'alt_title' : 'Filtering of Input Hierarchy',
            'contents' : [],
        } 

        txt = """
        > Filtering of hierarchy definitions (input vs optimised)
        Any group that is assigned no disorder in the optimised model is removed from the hierarchy definitions. 
        These new eff files are written to the output folder and can be used for further runs and may speed model optimisation. 
        """
        output['contents'] += self.format_summary(txt, classes=["square-corners-top"])

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
        output['contents'] += self.format_summary(txt)

        #        
        # Input Partitions
        #

        output_block = {
            'contents' : [],
        }
        output['contents'].append(output_block)

        output_block

        tab_set = {'type':'tabs', 'contents':[]}
        output_block['contents'].append(tab_set)

        # Add tab for each chain
        input_partitions = of.get('input_partitions_png',{})
        output_partitions = of.get('output_partitions_png',{})
        chain_ids = sorted(set(input_partitions.keys() + output_partitions.keys()))

        for c in chain_ids:

            p_in = input_partitions.get(c)
            p_out = output_partitions.get(c)
            tab = {
                'alt_title': 'Chain {}'.format(c),
                'contents' : [
                    {
                        'width':6,
                        'contents' : [{
                            'title':'Input Hierarchy',
                            'image':self.image(p_in),
                            }],
                        },
                    {
                        'width':6,
                        'contents' : [{
                            'title':'Filtered Hierarchy',
                            'image':self.image(p_out),
                            }],
                        },
                    ],
                }
            tab_set['contents'].append(tab)
        if tab_set['contents']:
            tab_set['contents'][0]['active'] = True

        return output