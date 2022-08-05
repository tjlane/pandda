import giant.logs as lg
logger = lg.getLogger(__name__)

import os
import pathlib as pl

from giant.html import (
    divs,
    formatters,
    )

from giant.html.summary import (
    HtmlSummary,
    )


class MakeMainPanddaHtmlPage(object):

    name = "Pandda Html Outputter"

    output_key = 'main_html'
    output_filename = 'pandda_main.html'
    template_name = 'pandda_main.html'

    def __init__(self,
        output_directory,
        dataset_html,
        analyse_html,
        inspect_html,
        ):

        self.output_directory = pl.Path(output_directory)
        self.output_path = (
            self.output_directory / self.output_filename
            )

        self.dataset_html = pl.Path(dataset_html)
        self.analyse_html = pl.Path(analyse_html)
        self.inspect_html = pl.Path(inspect_html)

    def __str__(self):

        s_ = (
            'Object: {name}\n'
            '| Output File:\n'
            '|\t{output_path}\n'
            '| Input Pages:\n'
            '|\t{dataset_html}\n'
            '|\t{analyse_html}\n'
            '|\t{inspect_html}\n'
            '`---->'
            ).format(
            name = self.name,
            output_path = str(self.output_path),
            dataset_html = str(self.dataset_html),
            analyse_html = str(self.analyse_html),
            inspect_html = str(self.inspect_html),
            )

        return s_.strip()

    def __call__(self):

        logger.subheading(
            'Writing output HTML: {}'.format(self.output_filename)
            )

        from pandda.html import HTML_ENV
        template = HTML_ENV.get_template(self.template_name)

        header_title = 'PanDDA Analysis Summary'
        body_header = self.get_body_header()

        # ===========================================================>
        # Construct the data object to populate the template
        output_data = {
            'header_title' : header_title,
            'body_header' : body_header,
            'pandda_dataset_html' : os.path.relpath( # use relative paths since within output file
                str(self.dataset_html.absolute()),
                str(self.output_path.absolute().parent),
                ),
            'pandda_analyse_html' : os.path.relpath(
                str(self.analyse_html.absolute()),
                str(self.output_path.absolute().parent),
                ),
            'pandda_inspect_html' : os.path.relpath(
                str(self.inspect_html.absolute()),
                str(self.output_path.absolute().parent),
                ),
            'contents' : [],
            }
        # ===========================================================>

        # Internal output objects
        output_data['contents'].extend(
            HtmlSummary.format_summary(
                """
                Output written to:
                <pre>{output_directory}</pre>
                Run output scripts, such as pandda.inspect from inside the output folder.
                """.format(
                    output_directory = str(self.output_directory.absolute()),
                )
            )
        )

        # Write out and format
        with open(str(self.output_path), 'w') as out_html:
            out_html.write(
                template.render(output_data)#.encode( "utf-8" )
                )

        logger(
            'Output HTML written to {}'.format(str(self.output_path))
            )

        return {
            self.output_key : str(self.output_path),
            'inspect_html' : str(self.inspect_html), # temporary hack?
            }

    def get_body_header(self):

        import pandda.resources

        format_citation = formatters.CitationFormatter(
            div_class = divs.Alert,
            width = 12,
            classes = ["bordered"],
            )

        d = {
            'title' : 'Pan-Dataset Density Analysis Results',
            'logo_image' : HtmlSummary.image(pandda.resources.SMALL_LOGO_PATH, embed=True),
            'logo_text' : pandda.resources.FORMAT_LOGO_TEXT.format('.analyse'),
            'citations' : [
                format_citation(
                    title = (
                        'A multi-crystal method for extracting obscured crystallographic '
                        'states from conventionally uninterpretable electron density.'
                        ),
                    journal = (
                        'Nature Communications, 8, 15123.'
                        ),
                    year = (
                        2017
                        ),
                    authors = (
                        'Pearce, N. M., Krojer, T., Bradley, A. R., Collins, P., '
                        'Nowak, R. P., Talon, R., Marsden, B. D., Kelm, S., '
                        'Shi, J., Deane, C. M., & von Delft, F.'
                        ),
                    link = (
                        'https://doi.org/10.1038/ncomms15123'
                        ),

                    ),
                format_citation(
                    title = (
                        'Proper modelling of ligand binding requires an ensemble of bound and unbound states.'
                        ),
                    journal = (
                        'Acta Cryst. Section D Structural Biology, 73(3), 256-266.'
                        ),
                    year = (
                        2017
                        ),
                    authors = (
                        'Pearce, N. M., Krojer, T., & von Delft, F.'
                        ),
                    link = (
                        'https://doi.org/10.1107/S2059798317003412'
                        ),
                    ),
                ],
            }

        return d
