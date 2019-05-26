import os, collections
import numpy, pandas
from libtbx import adopt_init_args
from pandemic.adp.html import HtmlSummary


class AnalyseResidualHtmlSummary(HtmlSummary):


    def __init__(self, task):
        adopt_init_args(self, locals())

    def main_summary(self):

        output = {
            'title' : 'Analysis of the difference between the target and the fitted ADPs/B-factors',
            'alt_title' : 'Residual Analysis',
            'contents' : [],
        }

        tab_set = {'type':'tabs', 'contents':[]}
        output['contents'].append(tab_set)

        tab_set['contents'].append(self.make_residuals())
        tab_set['contents'].append(self.make_correlations())

        if tab_set['contents']:
            tab_set['contents'][0]['active'] = True
        
        return [output]

    def make_residuals(self):

        of = self.task.result.output_files

        output = {
            'title' : 'Size of fitting residual [U(model) - U(target)] by chain',
            'alt_title' : 'Fitting Residual',
            'contents' : [],
        } 

        #        
        # Chain tab set
        #

        tab_set = {'type':'tabs', 'width':6, 'contents':[]}
        output['contents'].append(tab_set)
        for c, p in of.get('residuals_by_residue',{}).iteritems():
            tab = {
                'alt_title': 'Chain {}'.format(c),
                'contents' : [{'image':self.image(p)}],
                }
            tab_set['contents'].append(tab)
        if tab_set['contents']:
            tab_set['contents'][0]['active'] = True

        #
        # Residue tab set
        #

        tab_set = {'type':'tabs', 'width':6, 'contents':[]}
        output['contents'].append(tab_set)
        for l, p in of.get('residuals_by_dataset',{}).iteritems():
            tab = {
                'alt_title': l,
                'contents' : [{'image':self.image(p)}],
                }
            tab_set['contents'].append(tab)
        if tab_set['contents']:
            tab_set['contents'][0]['active'] = True

        #
        # Overall 
        #

        output['contents'].append({
            'width' : 6,
            #'title' : 'Fitting Residual vs Target B-factor',
            'image' : self.image(of.get('residuals_vs_bfactor')),
            })

        return output

    def make_correlations(self):

        of = self.task.result.output_files

        output = {
            'title' : 'Correlations with residual [U(model) - U(target)] by chain',
            'alt_title' : 'Fitting Residual Correlations',
            'contents' : [],
        } 

        tab_set = {'type':'tabs', 'contents':[]}
        output['contents'].append(tab_set)

        tab_hash = collections.OrderedDict()

        for l, l_dict in of['residual_correlations'].iteritems():

            for c, p in l_dict.iteritems():

                if c not in tab_hash:
                    tab_hash[c] = {
                        'title' : 'Correlations with residuals by level for chain {}'.format(c),
                        'alt_title' : 'Chain {}'.format(c),
                        'contents' : [{'type':'tabs', 'contents':[]}],
                        }

                tab_set2_contents = tab_hash[c]['contents'][0]['contents']

                tab_set2_contents.append({
                    #'title' : ' '.join([s.capitalize() for s in l.split(' ')]),
                    'alt_title': ' '.join([s.capitalize() for s in l.split(' ')]),
                    'contents' : [{'image' : self.image(p), 'width':8}],
                })
            
                tab_set2_contents[0]['active'] = True

        tab_set['contents'].extend(tab_hash.values())

        if tab_set['contents']:
            tab_set['contents'][0]['active'] = True

        return output

