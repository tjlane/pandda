import os, collections
import numpy, pandas

from pandemic.adp.html import (
    HtmlSummary,
    divs,
    )


class AnalyseResidualHtmlSummary(HtmlSummary):


    def __init__(self, task):
        self.task = task

    def main_summary(self):

        output = divs.Tab(alt_title='Residual Analysis')

        txt = """
        > Analysis of the difference between the target and the fitted ADPs/B-factors
        The fitting residual is the difference between target and fitted B-factors. The below tabs show the absolute size of this difference, as well as the correlation of the residuals to the different levels.
        """
        output.extend(self.format_summary(txt, classes=["square-corners-top"]))

        tab_set = output.append(divs.TabSet())
        tab_set.append(self.make_residuals())
        tab_set.append(self.make_correlations())
        tab_set.set_active()

        return [output]

    def make_residuals(self):

        of = self.task.result.output_files

        output = divs.Tab(alt_title='Fitting Residual Size')

        txt = """
        > Magnitude of the residuals.
        The size of the fitting residual is calculated as the RMS (root-mean-squared) value over all atoms.
        Large RMS values indicate that the fitted model fits poorly to the target B-factors.
        Graphs show the rms values by atom, but also grouped by dataset.
        """
        output.extend(self.format_summary(txt, classes=["square-corners-top"]))

        #
        # Chain tab set
        #

        chain_block = output.append(divs.Block())

        txt = """
        > Fitting residual by chain
        Large values for a particular atom may be due to a poorly resolved/modelled atom that has a non-physical B-factor.
        It can be expected that atoms at the end of sidechains (e.g. Lysine) will have larger RMSD values.
        """
        txt_block = divs.Block(
            width = 4,
            contents = self.format_summary(txt),
        )
        chain_block.append(txt_block)

        tab_set = chain_block.append(divs.TabSet(width=8))

        for c, p in of.get('residuals_by_residue',{}).items():
            tab = divs.Tab(
                alt_title = 'Chain {}'.format(c),
                contents = [divs.Block(image=self.image(p))],
            )
            tab_set.append(tab)
        tab_set.set_active()

        #
        # Dataset tab set
        #

        dataset_block = output.append(divs.Block())

        txt = """
        > Fitting residuals by dataset
        Large values for a particular dataset may be due to low dataset weights during optimisaton, e.g. for low resolution datasets.
        """
        txt_block = divs.Block(
            width = 4,
            contents = self.format_summary(txt),
        )
        dataset_block.append(txt_block)

        tab_set = dataset_block.append(divs.TabSet(width=8))

        for l, p in of.get('residuals_by_dataset',{}).items():
            tab = divs.Tab(
                alt_title = l,
                contents = [divs.Block(image=self.image(p))],
            )
            tab_set.append(tab)
        tab_set.set_active()

        #
        # Atom tab set
        #

        atom_type_block = output.append(divs.Block())

        txt = """
        > Fitting residuals by atom type
        Different atom types are better defined in the density that others, on average, and thus likely to give more physically-meaningful ADPs and therefore be fitted better.
        Due to this, comparison of each atom to others of the same type may more meaningfully indicate interesting outliers.
        """
        txt_block = divs.Block(
            width = 4,
            contents = self.format_summary(txt),
        )
        atom_type_block.append(txt_block)

        tab_set = divs.TabSet(width=8)
        atom_type_block.append(tab_set)

        for l, p in of.get('residuals_by_atom_type',{}).items():
            tab = divs.Tab(
                alt_title = l,
                contents = [divs.Block(image=self.image(p))],
            )
            tab_set.append(tab)
        tab_set.set_active()

        #
        # Overall
        #

        b_factor_block = output.append(divs.Block())

        txt = """
        > Fitting Residual vs Target B-factor
        The size of the residual plotted against the target B-factor of the input atoms.
        Generally, areas with larger B-factors are not significantly different to areas with smaller B-factors.
        """
        txt_block = divs.Block(
            width = 4,
            contents = self.format_summary(txt),
        )
        b_factor_block.append(txt_block)

        img_block = divs.Block(
            width = 8,
            image = self.image(of.get('residuals_vs_bfactor')),
        )
        b_factor_block.append(img_block)

        return output

    def make_correlations(self):

        of = self.task.result.output_files

        output = divs.Tab(alt_title = 'Fitting Residual Correlations')

        txt = """
        > Correlations between residual and model
        Large positive correlations between the model and the residual imply that the optimisation did not converge, or that the input model is incomplete.
        Random values across atoms indicate that the model had converged.
        """
        output.extend(self.format_summary(txt, classes=["square-corners-top"]))

        tab_set = output.append(divs.TabSet())

        tab_hash = collections.OrderedDict()

        for l, l_dict in of['residual_correlations'].items():

            for c, p in l_dict.items():

                # Changing the order of levels -> chains to chains -> levels
                if c not in tab_hash:
                    # Create new tab for this chain
                    c_tab = divs.Tab(alt_title = 'Chain {}'.format(c))
                    tab_hash[c] = c_tab
                    # Add text at top of tab
                    txt = """
                    > Correlations with residuals by level for chain {}
                    """.format(c)
                    c_tab.extend(self.format_summary(txt, ['square-corners-top']))

                # Extract for this chain
                c_tab = tab_hash[c]

                # Create block for this level
                l_block = c_tab.append(divs.Block())

                txt = """
                > Correlations for {} level
                """.format(l.title())
                l_block.extend(self.format_summary(txt, width=4))
                l_block.append(divs.Block(width=8, image=self.image(p)))

        tab_set.extend(tab_hash.values())
        tab_set.set_active()

        return output

