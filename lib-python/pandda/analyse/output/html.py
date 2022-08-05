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


class MakePanddaResultsHtml(object):

    output_key = 'analyse_html'
    output_filename = 'pandda_analyse.html'
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
        event_dicts,
        shell_dicts,
        dataset_dicts,
        output_files,
        ):

        contents = self.get_contents(
            event_dicts = event_dicts,
            shell_dicts = shell_dicts,
            dataset_dicts = dataset_dicts,
            output_files = output_files,
            )

        self.make_output(
            contents = contents,
            )

        return {self.output_key : str(self.output_path)}

    def get_contents(self, 
        event_dicts,
        shell_dicts,
        dataset_dicts,
        output_files,
        ):

        contents = []

        contents.extend(
            self.get_header()
            )

        main_block = divs.TabSet()
        contents.append(main_block)

        main_block.append(
            self.make_events_tab(
                event_dicts = event_dicts,
                test_dataset_keys = sorted([
                    k for k, v in dataset_dicts['test'].items()
                    if (v is True)
                    ]),
                output_graphs = output_files.get('graphs',{}),
                output_tables = output_files.get('tables',{}),
                )
            )
        main_block.set_active()

        main_block.append(
            self.make_shells_tab(
                shell_dicts = shell_dicts,
                output_graphs = output_files.get('graphs',{}),
                )
            )

        main_block.append(
            self.make_datasets_tab(
                dataset_dicts = dataset_dicts,
                all_dataset_keys = sorted(
                    dataset_dicts['high_resolution'].keys()
                    ),
                output_graphs = output_files.get('graphs',{}),
                output_tables = output_files.get('tables',{}),
                )
            )

        return contents

    def get_header(self):

        contents = []

        contents.append(
            divs.Block(
                title = "PanDDA results summary",
                fancy_title = True,
                )
            )

        return contents

    def make_events_tab(self,
        event_dicts,
        test_dataset_keys,
        output_graphs,
        output_tables,
        ):

        main_tab = divs.Tab(
            title = "Identified Events",
            alt_title = "Events",
            )

        main_tab.extend(
            self.get_events_header(
                event_dicts = event_dicts,
                test_dataset_keys = test_dataset_keys,
                )
            )

        main_tab.extend(
            self.get_events_main(
                event_dicts = event_dicts,
                output_graphs = output_graphs,
                output_tables = output_tables,
                )
            )

        return main_tab

    def get_events_header(self,
        event_dicts,
        test_dataset_keys,
        ):

        n_test_datasets = len(test_dataset_keys)
        n_with_events = len(set([e['dtag'] for e in event_dicts]))
        n_without_events = n_test_datasets - n_with_events

        block = divs.Block(
            classes = ['bordered'],
            contents = [
                objects.ProgressBar(
                    title = "Identified events in test datasets",
                    title_size = 4,
                    data = [
                        {
                            'label' : 'Tested datasets without events',
                            'value' : n_without_events,
                            'colour' : 'info',
                            },
                        {
                            'label' : 'Tested datasets with events',
                            'value' : n_with_events,
                            'colour' : 'success',
                            },
                        ],
                    add_counts = True,
                    add_percentages = True,
                    ),
                divs.Alert(
                    text = "Datasets without events: {}".format(n_without_events),
                    colour = "info",
                    width = 4,
                    ),
                divs.Alert(
                    text = "Datasets with events: {}".format(n_with_events),
                    colour = "success",
                    width = 4,
                    ),
                ],
            )

        return [block]

    def get_events_main(self,
        event_dicts,
        output_tables,
        output_graphs,
        ):

        site_counts = {}
        for e in event_dicts:
            s_num = e['site_num']
            site_counts.setdefault(s_num,0) 
            site_counts[s_num] += 1

        main_block = divs.Block(
            title = "Events Summary",
            classes = ['bordered'],
            contents = [
                divs.Block(
                    width = 6,
                    image = self.image(output_graphs.get('events_front')),
                    ),
                divs.Block(
                    width = 6,
                    image = self.image(output_graphs.get('events_back')),
                    ),
                divs.Block(
                    contents = [
                        objects.ProgressBar(
                            title = "Event distribution over sites",
                            data = [
                                {
                                    'label' : 'Site {}'.format(i),
                                    'value' : v,
                                    'colour' : ['default','info'][i%2],
                                    }
                                    for i,v in sorted(site_counts.items())
                                ],
                            ),
                        ],
                    ),
                divs.ScrollX(
                    contents = [
                        divs.Alert(
                            text = "Site {}: {} event(s)".format(i,v),
                            colour = "success",
                            width = 4,
                            )
                            for i,v in sorted(site_counts.items())
                        ],
                    ),            
                divs.ScrollX(
                    contents = [
                        divs.Block(
                            title = 'Site {}'.format(site_num),
                            text = 'Number of events: {}'.format(site_counts[site_num]),
                            image = self.image(site_image),
                            width = 5,
                            )
                        for site_num, site_image 
                        in sorted(output_graphs.get('site_events',{}).items())
                        ],
                    ),
                divs.Block(
                    contents = [
                        self.make_simple_block(
                            title = 'Estimated occupancies of events', 
                            image = output_graphs.get('event_fractions', None),
                            width = 6,
                            ),
                        self.make_simple_block(
                            title = 'Resolutions of identified events',
                            image = output_graphs.get('event_resolutions', None),
                            width = 6,
                            ),
                        ],
                    ),
                self.make_table_block(
                    title = 'Events Table', 
                    csv_file = output_tables.get('events_table', None), 
                    width = 12,
                    ),
                self.make_table_block(
                    title = 'Sites Table', 
                    csv_file = output_tables.get('sites_table', None), 
                    width = 12,
                    ),
                ],
            )

        return [main_block]

    def make_shells_tab(self,
        shell_dicts,
        output_graphs,
        ):

        main_tab = divs.Tab(
            title = "Map Analysis Shells",
            alt_title = "Shells",
            fancy_title = True,
            )

        main_tab.extend(
            self.get_shells_header(
                shell_dicts = shell_dicts,
                output_graphs = output_graphs,
                )
            )

        main_tab.extend(
            self.get_shells_main(
                shell_dicts = shell_dicts,
                output_graphs = output_graphs,
                )
            )

        return main_tab

    def get_shells_header(self,
        shell_dicts,
        output_graphs,
        ):

        n_shells = len(shell_dicts)
        n_datasets = sum([len(s['test']) for s in shell_dicts])

        min_res = min([s['resolution_high'] for s in shell_dicts])
        max_res = min([s['resolution_low'] for s in shell_dicts])

        main_block = divs.Block(
            contents = [
                divs.Alert(
                    text = (
                        "{n_datasets} datasets analysed in {n_shells} "
                        "shells from {min_res}-{max_res}A."
                        ).format(
                        n_datasets = n_datasets,
                        n_shells = n_shells,
                        min_res = min_res,
                        max_res = max_res,
                        ),
                    ),
                divs.Block(
                    contents = [
                        self.make_simple_block(
                            title = 'Map uncertainties of all datasets', 
                            image = output_graphs.get('map_uncertainties', None),
                            width = 6,
                            ),
                        self.make_simple_block(
                            title = 'Analysed Resolutions of all datasets',
                            image = output_graphs.get('analysed_resolution', None),
                            width = 6,
                            ),
                        ],
                    ),
                ],
            )

        return [main_block]

    def get_shells_main(self,
        shell_dicts,
        output_graphs,
        ):

        main_block = divs.TabSet(
            title = "Processing Shells",
            contents = [
                divs.Tab(
                    title = "Shell {shell_num}: {limit_high} - {limit_low}A".format(
                        shell_num = i_sd+1, 
                        limit_high = sd['resolution_high'],
                        limit_low = sd['resolution_low'],
                        ),
                    alt_title = sd['label'],
                    contents = [
                        self.get_shell_info_block(
                            info = sd,
                            ),
                        self.get_shell_images_block(
                            output_graphs = (
                                output_graphs.get('shells',{}).get(sd['label'],{})
                                ),
                            ),
                        ],
                    )
                    for i_sd, sd in enumerate(shell_dicts)
                ],
            )
        main_block.set_active()

        return [main_block]

    def get_shell_info_block(self, info):

        block = divs.Block(
            title = "Shell Info",
            contents = HtmlSummary.format_summary(
                (
                    "Resolutions: {limit_high} - {limit_low}\n"
                    "Map Resolution: {map_resolution}\n"
                    "Test Datasets: \n{test_datasets}\n"
                    "Training Datasets: \n{train_datasets}\n"
                    ).format(
                    limit_high = info['resolution_high'],
                    limit_low = info['resolution_low'],
                    map_resolution = info['map_resolution'],
                    test_datasets = ' '.join([
                        d_label.format(k=k) for k in info['test']
                        ]),
                    train_datasets = ' '.join([
                        d_label.format(k=k) for k in info['train']
                        ]),
                    )
                )
            )

        return block

    def get_shell_images_block(self, output_graphs):
        
        block = divs.Block(
            contents = [
                self.make_pretty_block(
                    title = "Shell Statistical Map Distributions",
                    image = output_graphs.get('statistical_map_distribution'),
                    text = (
                        "Distribution of map values for the mean (mu) and variation (sigma-a) "
                        "maps. "
                        ),
                    ),
                self.make_pretty_block(
                    title = "Shell Dataset Uncertainties",
                    image = output_graphs.get('uncertainty_distribution'),
                    text = (
                        "Uncertainties of training and test datasets."
                        ),
                    ),
                ],
            )

        return block

    def make_datasets_tab(self,
        dataset_dicts,
        all_dataset_keys,
        output_graphs,
        output_tables,
        ):

        main_tab = divs.Tab(
            title = "Dataset-by-Dataset Summary",
            alt_title = "Datasets",
            )

        main_tab.extend(
            self.get_datasets_header(
                dataset_dicts = dataset_dicts,
                )
            )

        main_tab.extend(
            self.get_datasets_main(
                dataset_dicts = dataset_dicts,
                all_dataset_keys = all_dataset_keys,
                output_graphs = output_graphs,
                output_tables = output_tables,
                )
            )

        return main_tab

    def get_datasets_header(self,
        dataset_dicts,
        ):

        train_keys = [
            k for k, v in dataset_dicts['train'].items()
            if (v is True)
            ]
        test_keys = [
            k for k, v in dataset_dicts['test'].items()
            if (v is True)
            ]
        
        n_test = len(test_keys)
        n_test_and_train = len(set(test_keys).intersection(train_keys))
        n_test_only = len(test_keys) - n_test_and_train
        n_train_only = len(train_keys) - n_test_and_train

        block = divs.Block(
            classes = ['bordered'],
            contents = [
                objects.ProgressBar(
                    title = "Test/Train Datasets",
                    title_size = 4,
                    text = (
                        "Numbers of datasets used for training the statistical model "
                        "(building a picture of the ground state of the crystal), "
                        "testing against the statistical model (looking for events), "
                        "or both."
                        ),
                    data = [
                        {
                            'label' : 'Training only',
                            'value' : n_train_only,
                            'colour' : 'info',
                            },
                        {
                            'label' : 'Training & testing',
                            'value' : n_test_and_train,
                            'colour' : 'warning',
                            },
                        {
                            'label' : 'Testing only',
                            'value' : n_test_only,
                            'colour' : 'success',
                            },
                        ],
                    add_counts = True,
                    add_percentages = True,
                    ),
                divs.Alert(
                    text = "Train-only datasets: {}".format(n_train_only),
                    colour = "info",
                    width = 4,
                    ),
                divs.Alert(
                    text = "Test&Train datasets: {}".format(n_test_and_train),
                    colour = "warning",
                    width = 4,
                    ),
                divs.Alert(
                    text = "Test-only datasets: {}".format(n_test),
                    colour = "success",
                    width = 4,
                    ),
                ],
            )

        return [block]

    def get_datasets_main(self,
        dataset_dicts,
        all_dataset_keys,
        output_graphs,
        output_tables,
        ):

        all_contents = []

        for dkey in all_dataset_keys:

            d_info = {
                k : d.get(dkey)
                for k, d in dataset_dicts.items()
                }

            d_contents = [
                self.get_dataset_info_block(
                    info = d_info,
                    ),
                self.get_dataset_base_images_block(
                    dataset_key = dkey,
                    output_graphs = output_graphs,
                    ),
                ]

            if d_info['train'] is True:
                pass

            if d_info['test'] is True:
                d_contents.append(
                    self.get_dataset_test_images_block(
                        dataset_key = dkey,
                        output_graphs = output_graphs,
                        )
                    )

            if output_graphs.get('bdc_estimation',{}).get(dkey):
                d_contents.append(
                    self.get_dataset_event_images_block(
                        dataset_key = dkey,
                        output_graphs = output_graphs,
                        )
                    )

            all_contents.append(
                divs.Tab(
                    title = dkey, 
                    alt_title = dkey,
                    contents = d_contents,
                    )
                )

        main_block = divs.TabSet(
            title = "Processed Datasets",
            contents = all_contents,
            )
        main_block.set_active()

        return [main_block]

    def get_dataset_info_block(self, info):

        block = divs.Block(
            title = "Dataset Info",
            contents = HtmlSummary.format_summary(
                text = (
                    "High Resolution: {high_res}\n"
                    "Low Resolution: {low_res}\n"
                    "R-free / R-work: {r_free} / {r_work}\n"
                    "Space Group: {space_group}\n"
                    "Unit cell: {unit_cell}\n"
                    "\n"
                    "Test Dataset: {test}\n"
                    "Training Dataset: {train}\n"
                    "Analysed Resolution: {analysed_resolution}\n"
                    "Map Uncertainty: {map_uncertainty}\n"
                    ).format(
                    high_res = info['high_resolution'],
                    low_res = info['low_resolution'],
                    r_free = info['r_free'],
                    r_work = info['r_work'],
                    space_group = info['space_group'],
                    unit_cell = info['unit_cell'],
                    test = info['test'],
                    train = info['train'],
                    analysed_resolution = info.get('analysed_resolution', 'n/a'),
                    map_uncertainty = info.get('map_uncertainty', 'n/a')
                    ),
                width = 6,
                ),
            )

        return block

    def get_dataset_base_images_block(self, dataset_key, output_graphs):

        block = divs.Block(
            contents = [
                self.make_pretty_block(
                    title = "Dataset Wilson Plot",
                    image = output_graphs.get('wilson_plot',{}).get(dataset_key),
                    text = (
                        "Wilson plot of dataset reflection data and the reference dataset."
                        ),
                    ),
                ],
            )

        return block

    def get_dataset_test_images_block(self, dataset_key, output_graphs):

        block = divs.Block(
            contents = [
                divs.Alert(
                    title = "Evaluations graphs (for test datasets)",
                    ),
                self.make_pretty_block(
                    title = "Dataset Z-Map Distribution",
                    image = output_graphs.get('map_distribution',{}).get(dataset_key),
                    text = (
                        "Distribution of Z-values in the dataset map. After "
                        "normalisation, the Z-values should be approximately "
                        "normally-distributed."
                        ),
                    ),
                self.make_pretty_block(
                    title = "Dataset QQ-plot",
                    image = output_graphs.get('qq_plot',{}).get(dataset_key),
                    text = (
                        "Quantile-quantile plot used to calculate the noise of the dataset."
                        ),
                    ),
                self.make_pretty_block(
                    title = "Map Scatter",
                    image = output_graphs.get('map_scatter',{}).get(dataset_key),
                    text = (
                        "Scatter plots of the dataset map against the mean map. "
                        "The deviation from the diagonal is indicative of the noise in "
                        "the dataset. The sorted maps should show an almost perfect straight "
                        "line -- any deviation from the diagonal here indicates a significant "
                        "problem with the dataset/the dataset map."
                        ),
                    ),
                ],
            )

        return block

    def get_dataset_event_images_block(self, dataset_key, output_graphs):

        block = divs.Block(
            contents = [
                divs.Alert(
                    title = "Identified Event Graphs",
                    ),
                divs.ScrollX(
                    title = "Event Background Estimation",
                    contents = [
                        divs.Block(
                            title = 'Event {}'.format(event_num),
                            image = self.image(event_image),
                            width = 5,
                            )
                        for event_num, event_image 
                        in sorted(output_graphs.get('bdc_estimation',{}).get(dataset_key).items())
                        ],
                    ),
                ],
            )

        return block

    def make_simple_block(self, title, image, text=None, width=12):

        block = divs.Block(
            width = width,
            title = title,
            image = self.image(image),
            contents = (
                HtmlSummary.format_summary(text, type="block")
                if (text is not None) else []
                )
            )

        return block

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
            '<th></th>','<th>Dataset</th>'
            ).replace(
            'border="1" ', ''
            )

        block = divs.Alert(
            title = title,
            width = width,
            table = table_html,
            )

        return block

    def make_pretty_block(self, title, image, text, width=12, text_width=6, image_width=6):

        block = divs.Block(
            width = width,
            title = title,
            contents = [
                divs.Block(
                    contents = HtmlSummary.format_summary(text),
                    width = text_width,
                    ),
                divs.Block(
                    image = self.image(image),
                    width = image_width,
                    ),
                ],
            )

        return block

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
