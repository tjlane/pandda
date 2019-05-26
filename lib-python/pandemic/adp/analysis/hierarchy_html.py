import os, collections
import numpy, pandas
from libtbx import adopt_init_args
from pandemic.adp.html import HtmlSummary


class AssessHierarchyGroupsHtmlSummary(HtmlSummary):


    def __init__(self, task):
        adopt_init_args(self, locals())

    def main_summary(self):

        output = {
            'title' : 'Analysis of hierarchy definitions (input vs optimised)',
            'alt_title' : 'Hierarchy Definition',
            'contents' : [],
        }

        tab_set = {'type':'tabs', 'contents':[]}
        output['contents'].append(tab_set)

        tab_set['contents'].append(self.make_panel())

        if tab_set['contents']:
            tab_set['contents'][0]['active'] = True
        
        return [output]

    def make_panel(self):

        of = self.task.result.output_files

        output = {
            'title' : 'XXX',
            'alt_title' : 'XXX',
            'contents' : [],
        } 

        #        
        # Output Partitions
        #

        tab_set = {'type':'tabs', 'width':6, 'contents':[]}
        output['contents'].append(tab_set)
        for c, p in of.get('level_partitions',{}).iteritems():
            tab = {
                'alt_title': 'Chain {}'.format(c),
                'contents' : [{'image':self.image(p)}],
                }
            tab_set['contents'].append(tab)
        if tab_set['contents']:
            tab_set['contents'][0]['active'] = True

        return output