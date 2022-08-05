import giant.logs as lg
logger = lg.getLogger(__name__)

import pathlib as pl
import numpy as np
import pandas as pd

from giant.html import (
    divs, 
    objects,
    formatters,
    )

from giant.html.summary import (
    ImageEmbedder,
    HtmlSummary,
    )

d_label = '<span class="label label-info">{k}</span>'


class MakePanddaDatasetSummaryHtml(object):

    output_filename = 'pandda_initial.html'
    template_name = 'pandda_page.html'

    def __init__(self,
        output_directory,
        ):

        self.output_directory = pl.Path(output_directory)
        self.output_path = (
            self.output_directory / self.output_filename
            )

        self.image = ImageEmbedder(
            embed = False,
            relative_to = str(self.output_path.parent),
            )
        
    def __call__(self,
        datasets,
        output_info,
        output_graphs,
        output_tables,
        ):

        contents = self.get_contents(
            datasets = datasets,
            output_info = output_info,
            output_graphs = output_graphs,
            output_tables = output_tables,
            )

        self.make_output(
            contents = contents,
            )

        return self.output_path

    def get_contents(self, 
        datasets,
        output_info,
        output_graphs,
        output_tables,
        ):

        contents = []

        contents.extend(
            self.get_header()
            )

        contents.extend(
            self.get_input_summary(
                info_dict = output_info,
                )
            )

        contents.extend(
            self.get_mcd_summary(
                datasets = datasets, 
                output_graphs = output_graphs,
                output_tables = output_tables,
                )
            )

        return contents

    def get_header(self):

        contents = []

        contents.append(
            divs.Block(
                title = "Input Dataset Summary",
                fancy_title = True,
                )
            )

        return contents

    def get_input_summary(self,
        info_dict,
        ):

        empty_directories = sorted(
            info_dict['empty_directories']
            )

        filtered_files = sorted(
            info_dict['filtered_files'].items()
            )

        rejected_datasets = sorted(
            info_dict['rejected_datasets'].items()
            )

        loaded_datasets = sorted(
            info_dict['loaded_datasets']
            )

        n_empty = len(empty_directories)
        n_filtered = len(filtered_files)
        n_rejected = len(rejected_datasets)
        n_loaded = len(loaded_datasets)

        ###

        panel = divs.Panel(title='Input Files/Folders/Datasets')

        panel.extend([
            objects.ProgressBar(
                title = "Input dataset summary",
                title_size = 4,
                data = [
                    {
                        'label' : 'Empty Directories',
                        'value' : n_empty,
                        'colour' : 'danger',
                        },
                    {
                        'label' : 'Filtered Files',
                        'value' : n_filtered,
                        'colour' : 'info',
                        },
                    {
                        'label' : 'Rejected Datasets',
                        'value' : n_rejected,
                        'colour' : 'warning',
                        },
                    {
                        'label' : 'Successfully Loaded',
                        'value' : n_loaded,
                        'colour' : 'success',
                        }
                    ],
                add_counts = True,
                add_percentages = True,
                ),
            divs.Alert(
                text = "Empty Directories: {}".format(n_empty),
                width = 3,
                colour = "danger",
                ),
            divs.Alert(
                text = "Filtered Files: {}".format(n_filtered),
                width = 3,
                colour = "info",
                ),
            divs.Alert(
                text = "Rejected Datasets: {}".format(n_rejected),
                width = 3,
                colour = "warning",
                ),
            divs.Alert(
                text = "Processed Datasets: {}".format(n_loaded),
                width = 3,
                colour = "success",
                ),
            ])

        # 

        if empty_directories: 
            panel.append(
                divs.Panel(
                    title = "Empty Directories",
                    show = False,
                    contents = [
                        divs.Block(
                            text = '<br>'.join(empty_directories),
                            ),
                        ],
                    )
                )

        if filtered_files: 

            filtered_files_counts = {}
            for dkey, (reason, files) in filtered_files:
                filtered_files_counts.setdefault(reason,[]).append(dkey)

            block = panel.append(
                divs.Panel(
                    title = "Ignored Datasets",
                    show = True,
                    )
                )

            for reason, keys in sorted(filtered_files_counts.items()):
                block.append(
                    divs.Block(
                        title = "{r} ({n} datasets)".format(
                            r=reason, n=len(keys),
                            ),
                        text = ' '.join([
                            d_label.format(k=k)
                            for k in keys
                            ]),
                        )
                    )

        if rejected_datasets:

            rejected_datasets_counts = {}
            for dkey, reason in rejected_datasets:
                rejected_datasets_counts.setdefault(reason,[]).append(dkey)

            block = panel.append(
                divs.Panel(
                    title = "Rejected Datasets",
                    show = True,
                    )
                )

            for reason, keys in sorted(rejected_datasets_counts.items()):
                block.append(
                    divs.Block(
                        title = "Rejection Reason: {r} ({n} datasets)".format(
                            r=reason, n=len(keys),
                            ),
                        text = ' '.join([
                            d_label.format(k=k)
                            for k in keys
                            ]),
                        )
                    )

        if loaded_datasets: 

            panel.append(
                divs.Panel(
                    title = "Loaded Datasets",
                    show = False,
                    contents = [
                        divs.Block(
                            text = ' '.join([
                                d_label.format(k=k)
                                for k in loaded_datasets
                                ]),
                            ),
                        ],
                    )
                )

        return [panel]

    def get_mcd_summary(self, 
        datasets, 
        output_graphs,
        output_tables,
        ):

        width = 6

        main_block = divs.Panel(
            title = "Loaded Datasets Summary",
            contents = [
                self.make_standard_block(
                    title = (
                        "Distribution of dataset resolutions."
                        ),
                    width = width,
                    image = output_graphs.get('resolution_distribution'),
                    ),
                self.make_standard_block(
                    title = (
                        "Distributions of dataset R-factors."
                        ),
                    width = width,
                    image = output_graphs.get('rfactor_distribution'),
                    ),
                self.make_standard_block(
                    title = (
                        "Distributions of dataset unit cells."
                        ),
                    width = width,
                    image = output_graphs.get('unit_cell_distribution'),
                    ),
                self.make_table_block(
                    title = (
                        "Dataset information"
                        ),
                    width = 12,
                    csv_file = output_tables.get('dataset_info'),
                    ),
                ],
            )
        
        return [main_block]

    def make_standard_block(self, title, image, text=None, width=12):

        main_block = divs.Block(
            width = width,
            title = title,
            image = self.image(image),
            contents = (
                HtmlSummary.format_summary(text, type="block")
                if (text is not None) else []
                )
            )

        return main_block

    def make_table_block(self, title, csv_file, width=12):

        table = pd.read_csv(
            csv_file,
            index_col=False,
            )

        table_html = table.to_html(
            index = False,
            bold_rows = False,
            na_rep = '',
            classes = ['table table-striped table-hover datatable nowrap'],
            ).replace(
            'border="1" ', ''
            )

        main_block = divs.Alert(
            title = title,
            width = width,
            table = table_html,
            )

        return main_block

    def make_output(self,
        contents,
        ):

        logger.subheading(
            'Writing output HTML: {}'.format(self.output_filename)
            )

        if not self.output_directory.exists():
            self.output_directory.mkdir(parents=True)

        from pandda.html import HTML_ENV
        template = HTML_ENV.get_template(self.template_name)

        header_title = 'Pandda Input Datasets Summary'

        # ===========================================================>
        # Construct the data object to populate the template
        output_data = {
            'header_title' : header_title,
            #'body_header' : body_header,
            'contents' : contents,
            }
        # Jsons
        #if json_plots:
        #    output_data['json_plots'] = json_plots
        # ===========================================================>

        # Write out and format
        with open(str(self.output_path), 'w') as out_html:
            out_html.write(
                template.render(output_data)#.encode( "utf-8" )
                )

        logger(
            'Output HTML written to {}'.format(str(self.output_path))
            )
