import giant.logs as lg
logger = lg.getLogger(__name__)

import pathlib as pl

from giant.html import (
    divs,
    formatters,
    )

from giant.html.summary import (
    HtmlSummary,
    HtmlSummaryCollator,
    HtmlSummaryConcatenator,
    as_html_summary_maybe,
    as_html_summaries_maybe,
    )


class WriteHtmlSummaryTask(object):

    output_filename = 'results.html'
    template_name = 'adp_summary.html'

    def __init__(self,
        output_directory,
        ):
        
        self.output_directory = pl.Path(output_directory)
        self.output_path = (self.output_directory / self.output_filename)

    def run(self,
        html_objects,
        ):

        # Things for the first page
        short_summaries = []
        main_summaries = []
        json_plots = []

        # Compile output lists
        for obj in html_objects:
            if hasattr(obj, 'short_summary'):
                short_summaries.extend(obj.short_summary())
            if hasattr(obj, 'main_summary'):
                main_summaries.extend(obj.main_summary())
            if hasattr(obj, 'json_plots'):
                json_plots.extend(obj.json_plots())

        self.make_output(
            overview_objects = short_summaries,
            tab_objects = main_summaries,
            json_plots = json_plots,
            )

        return self

    def make_output(self,
        overview_objects,
        tab_objects,
        json_plots = [],
        ):

        logger.subheading(
            'Writing output HTML: {}'.format(self.output_filename)
            )

        from pandemic.html import HTML_ENV
        template = HTML_ENV.get_template(self.template_name)

        header_title = 'PanDEMIC ADP Summary'
        body_header = self.get_body_header()

        # ===========================================================>
        # Construct the data object to populate the template
        output_data = {
            'header_title' : header_title,
            'body_header' : body_header,
            'contents' : [],
            }
        # Jsons
        if json_plots:
            output_data['json_plots'] = json_plots
        # ===========================================================>

        # Internal output objects
        output_data['contents'].extend(
            HtmlSummary.format_summary(
                """
                Output written to:
                <pre>{output_directory}</pre>
                If the output folder has been moved, scripts (e.g. pymol scripts) must be run from inside the folder containing the script.
                """.format(
                    output_directory = str(self.output_directory.absolute()),
                )
            )
        )

        # Create overview tab
        tab_set = divs.TabSet()
        tab_set.append(self.create_overview_tab(objects=overview_objects))
        tab_set.set_active()

        # Create other tabs
        tab_set.extend(tab_objects)

        output_data['contents'].append(tab_set)

        # Write out and format
        with open(str(self.output_path), 'w') as out_html:
            out_html.write(
                template.render(output_data)#.encode( "utf-8" )
                )

        logger(
            'Output HTML written to {}'.format(str(self.output_path))
            )

    def get_body_header(self):

        import pandemic.resources

        format_citation = formatters.CitationFormatter()

        d = {
            'title' : 'Hierarchical Disorder Characterisation Summary',
            'logo_image' : HtmlSummary.image(pandemic.resources.SMALL_LOGO_PATH, embed=True),
            'logo_text' : pandemic.resources.FORMAT_LOGO_TEXT.format('.adp'),
            'citations' : [
                format_citation(
                    title = 'A method for the intuitive extraction of macromolecular dynamics from structural disorder.',
                    journal = 'In submission',
                    authors = ['N. Pearce', 'P. Gros'],
                    year = "2021",
                    link = None,
                    ),                
                ],
            }

        return d

    def create_overview_tab(self, objects):

        tab = divs.Tab(
            id = 'overview',
            title = 'Hierarchical Disorder Parameterisation Overview',
            alt_title = 'Overview',
            contents = objects,
        )

        return tab

