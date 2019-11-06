import os, collections
import numpy, pandas
from libtbx import adopt_init_args
from pandemic.adp.html import HtmlSummary


class AnalyseResidualHtmlSummary(HtmlSummary):


    def __init__(self, task):
        adopt_init_args(self, locals())

    def main_summary(self):

        output = {
            'alt_title' : 'Residual Analysis',
            'contents' : [],
        }

        txt = """> Analysis of the difference between the target and the fitted ADPs/B-factors
        The fitting residual is the difference between target and fitted B-factors. The below tabs show the absolute size of this difference, as well as the correlation of the residuals to the different levels.
        """
        output['contents'].extend(self.format_summary(txt))

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
            'title' : '',
            'alt_title' : 'Fitting Residual Size',
            'contents' : [],
        } 

        txt = """The size of the fitting residual is calculated as the RMS (root-mean-squared) value over all atoms.
        Large RMS values indicate that the fitted model fits poorly to the target B-factors.
        Graphs show the rms values by atom, but also grouped by dataset.
        """
        output['contents'].extend(self.format_summary(txt))

        #        
        # Chain tab set
        #

        chain_block = {'width':12, 'contents' : []}
        output['contents'].append(chain_block)

        txt = """
        > Fitting residual by chain
        Large values for a particular atom may be due to a poorly resolved/modelled atom that has a non-physical B-factor. It can be expected that atoms at the end of sidechains (e.g. Lysine) will have larger RMSD values.
        """
        txt_block = {'width':4, 'contents' : self.format_summary(txt)}
        chain_block['contents'].append(txt_block)

        tab_set = {
            'type':'tabs', 
            'width' : 8,
            'contents':[],
            }
        chain_block['contents'].append(tab_set)

        for c, p in of.get('residuals_by_residue',{}).iteritems():
            tab = {
                'alt_title': 'Chain {}'.format(c),
                'contents' : [{'image':self.image(p)}],
                }
            tab_set['contents'].append(tab)
        if tab_set['contents']:
            tab_set['contents'][0]['active'] = True

        #
        # Dataset tab set
        #

        dataset_block = {'width':12, 'contents' : []}
        output['contents'].append(dataset_block)

        txt = """
        > Fitting residuals by dataset
        Large values for a particular dataset may be due to low dataset weights during optimisaton, e.g. for low resolution datasets.
        """
        txt_block = {'width':4, 'contents' : self.format_summary(txt)}
        dataset_block['contents'].append(txt_block)

        tab_set = {
            'type' : 'tabs', 
            'width' : 8, 
            'contents':[],
            }
        dataset_block['contents'].append(tab_set)

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

        b_factor_block = {'width':12, 'contents' : []}
        output['contents'].append(b_factor_block)

        txt = """> Fitting Residual vs Target B-factor
        """
        txt_block = {'width':4, 'contents' : self.format_summary(txt)}
        b_factor_block['contents'].append(txt_block)

        img_block = {
            'width' : 8, 
            'image' : self.image(of.get('residuals_vs_bfactor')),
            }
        b_factor_block['contents'].append(img_block)

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

