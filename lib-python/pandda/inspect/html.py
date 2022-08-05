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
    )

from .exceptions import (
    MissingFile
    )

d_label = '<span class="label label-info">{k}</span>'


class MakePanddaInspectHtml(object):

    output_filename = 'pandda_inspect.html'
    template_name = 'pandda_page.html'

    def __init__(self, 
        output_directory,
        ):

        self.output_directory = pl.Path(output_directory)

        self.input_files = self.build_input_files_dict()

        self.input_path_prefix = self.find_input_prefix(
            filepath = self.input_files['html']['main_html'],
            )

        self.output_path = pl.Path(
            self.get_path(['html','inspect_html'])
            )

        self.image = ImageEmbedder(
            embed = False,
            relative_to = str(self.output_path.parent),
            )

    def __call__(self,
        inspector,
        output_files,
        ):

        results = self.unpack_inspector(
            inspector = inspector,
            )

        contents = self.get_contents(
            results = results,
            output_files = output_files,
            )

        self.make_output(
            contents = contents,
            )

        return self.output_path

    def build_input_files_dict(self):

        import json

        json_file = (
            self.output_directory / 'results.json'
            )

        if not json_file.exists(): 
            raise MissingFile('Json results file not found: {}'.format(str(json_file)))

        json_string = open(str(json_file), 'r').read()

        pandda_log_obj = json.loads(json_string)

        input_files = pandda_log_obj['output_files']

        return input_files

    def find_input_prefix(self, filepath):

        filepath = pl.Path(filepath)

        prefix = filepath.parent

        for _ in filepath.parts:

            test_f = filepath.relative_to(prefix)

            if (self.output_directory / test_f).exists():
                break

            prefix = prefix.parent
        
        return prefix

    def get_path(self, input_keys):

        assert isinstance(input_keys, list)

        d = self.input_files
        for k in input_keys:
            d = d[k]

        assert not hasattr(d, 'keys')

        output_path = (
            self.output_directory / pl.Path(d).relative_to(
                self.input_path_prefix
                )
            )

        return str(output_path)

    def unpack_inspector(self,
        inspector,
        ):

        results = {
            'event_table' : inspector.tables.events.table.reset_index(),
            'site_table' : inspector.tables.sites.table.reset_index(),
        }

        return results

    def get_contents(self, 
        results,
        output_files,
        ):

        contents = []

        contents.extend(
            self.get_header()
            )

        contents.extend(
            self.get_main(
                results = results,
                output_files = output_files,
                )
            )

        return contents

    def get_header(self):

        contents = [
            divs.Block(
                title = "PanDDA inspection summary",
                fancy_title = True,
                ),
            ]

        return contents

    def get_main(self,
        results,
        output_files,
        ):

        event_table = results['event_table']

        event_table_dict = self.build_event_table_summary_dict(
            event_table = event_table,
            )

        ###

        contents = []

        contents.extend(
            self.get_inspection_progress_bars(
                summary_dict = event_table_dict,
                )
            )

        contents.extend(
            self.get_inspection_site_summary(
                summary_dict = event_table_dict,
                output_files = output_files,
                )
            )

        contents.extend(
            self.get_inspection_event_summary(
                event_table = event_table,
                )
            )

        return contents

    def get_inspection_progress_bars(self,
        summary_dict,
        ):

        block1 = divs.Block(
            classes = ['bordered'],
            contents = [
                objects.ProgressBar(
                    title = "Inspection Progress",
                    title_size = 4,
                    data = [
                        {
                            'label' : 'Fitted',
                            'value' : summary_dict['n_fitted'],
                            'colour' : 'success',
                            },
                        {
                            'label' : 'Unviewed',
                            'value' : summary_dict['n_unviewed'],
                            'colour' : 'info',
                            },
                        {
                            'label' : 'No Ligand Fitted',
                            'value' : summary_dict['n_empty'],
                            'colour' : 'danger',
                            },
                        ],
                    add_counts = True,
                    add_percentages = False,
                    ),
                divs.Alert(
                    text = "Total events: {}".format(
                        summary_dict['n_blobs']
                        ),
                    colour = "info",
                    width = 12,
                    ),
                divs.Alert(
                    text = "Fitted ligands: {}".format(
                        summary_dict['n_fitted']
                        ),
                    colour = "success",
                    width = 4,
                    ),
                divs.Alert(
                    text = "Unviewed: {}".format(
                        summary_dict['n_unviewed']
                        ),
                    colour = "info",
                    width = 4,
                    ),
                divs.Alert(
                    text = "No ligand fitted: {}".format(
                        summary_dict['n_empty']
                        ),
                    colour = "danger",
                    width = 4,
                    ),
                ],
            )

        block2 = divs.Block(
            classes = ['bordered'],
            contents = [
                divs.Alert(
                    text = "Datasets with ligands: {}".format(
                        summary_dict['n_datasets_w_hit']
                        ),
                    colour = "info",
                    width = 6,
                    ),
                divs.Alert(
                    text = "Sites with ligands: {}".format(
                        summary_dict['n_sites_w_hit']
                        ),
                    colour = "info",
                    width = 6,
                    ),
                ],
            )

        block3 = divs.Block(
            classes = ['bordered'],
            contents = (
                [
                    objects.ProgressBar(
                        title = "Fitted ligands",
                        title_size = 4,
                        data = [
                            {
                                'label' : 'High confidence',
                                'value' : summary_dict['n_high_confidence'],
                                'colour' : 'success',
                                },
                            {
                                'label' : 'Medium confidence',
                                'value' : summary_dict['n_medium_confidence'],
                                'colour' : 'info',
                                },
                            {
                                'label' : 'Low confidence',
                                'value' : summary_dict['n_low_confidence'],
                                'colour' : 'danger',
                                },
                            {
                                'label' : 'Unranked',
                                'value' : summary_dict['n_unranked'],
                                'colour' : 'info',
                                },
                            ],
                        add_counts = True,
                        add_percentages = False,
                        ),
                    divs.Alert(
                        text = "Total ligands: {}".format(
                            summary_dict['n_fitted']
                            ),
                        colour = "info",
                        width = 12,
                        ),
                    divs.Alert(
                        text = "High Confidence: {}".format(
                            summary_dict['n_high_confidence']
                            ),
                        colour = "success",
                        width = 4,
                        ),
                    divs.Alert(
                        text = "Medium Confidence: {}".format(
                            summary_dict['n_medium_confidence']
                            ),
                        colour = "warning",
                        width = 4,
                        ),
                    divs.Alert(
                        text = "Low Confidence: {}".format(
                            summary_dict['n_low_confidence']
                            ),
                        colour = "danger",
                        width = 4,
                        ),
                    divs.Alert(
                        text = "Unfitted, marked interesting: {}".format(
                            summary_dict['n_interesting_unfitted']
                            ),
                        colour = "danger",
                        width = 6,
                        ),
                    ]
                if (summary_dict['n_fitted'] > 0) else 
                [
                    divs.Alert(
                        text = "Nothing to show here yet: no fitted ligands"
                        )
                    ]
                ),
            )

        return [block1, block2, block3]

    def get_inspection_site_summary(self,
        summary_dict,
        output_files,
        ):

        site_counts = summary_dict['n_hits_per_site']
        site_images = output_files['site_events']

        blocks = [
            divs.Block(
                width = 6,
                text = 'Sites and Events (front)',
                image = self.get_path(['graphs','events_front']),
                ),
            divs.Block(
                width = 6,
                text = 'Sites and Events (back)',
                image = self.get_path(['graphs','events_back']),
                ),
            (
                objects.ProgressBar(
                    title = "Fitted ligands distribution over sites",
                    data = [
                        {
                            'label' : 'Site {}'.format(i),
                            'value' : v,
                            'colour' : ['default','info'][i%2],
                            }
                            for i,v in sorted(site_counts.items())
                        ],
                    add_counts = True,
                    add_percentages = True,
                    )
                if sum(site_counts.values()) > 0 else
                divs.Alert(
                    text = 'Nothing to show here yet: no fitted ligands',
                    )
                ),
            divs.ScrollX(
                contents = [
                    divs.Block(
                        title = 'Site {}'.format(site_num),
                        text = 'Number of events: {}'.format(n_events),
                        image = self.image(site_images[site_num]),
                        width = 5,
                        )
                    for site_num, n_events
                    in sorted(site_counts.items())
                    ],
                ),
            ]

        return blocks

    def get_inspection_event_summary(self,
        event_table,
        ):

        column_labels = [] # subsample the columns

        # ['Dataset','Viewed','Interesting','Lig. Placed','Event','Site','1 - BDC','Z-Peak','Map Res.','Map Unc.','Confidence','Comment','']

        blocks = [
            self.make_table_block(
                title = 'Events/Fitted Ligands Table', 
                table = event_table, 
                width = 12,
                ),
            ]

        return blocks

    def build_event_table_summary_dict(self,
        event_table,
        ):

        # Input counts

        n_blobs = len(
            event_table.index
            )

        n_sites = len(
            set(event_table['site_num'])
            )

        n_datasets = len(
            # this assumes the index is (dtag, e_idx)
            set(event_table['dtag'])
            )

        # Blobs Inspected/Modelled/Empty

        n_fitted = sum(
            event_table['Ligand Placed']
            )

        n_viewed = sum(
            event_table['Viewed']
            )

        n_empty = (
            n_viewed - n_fitted
            )

        n_unviewed = (
            n_blobs - n_viewed
            )

        # Interesting unfitted (bookmarked)

        n_interesting_unfitted = sum(
            event_table["Interesting"][event_table['Ligand Placed']==False]
            )

        # Confidence of models

        n_high_confidence = sum(
            event_table["Ligand Placed"][event_table['Ligand Confidence']=='High']
            )

        n_medium_confidence = sum(
            event_table["Ligand Placed"][event_table['Ligand Confidence']=='Medium']
            )
        
        n_low_confidence = sum(
            event_table["Ligand Placed"][event_table['Ligand Confidence']=='Low']
            )

        n_unranked = (
            n_fitted - n_high_confidence - n_medium_confidence - n_low_confidence
            )

        # Datasets/sites with hits

        try:    
            n_datasets_w_hit = len(
                set(zip(*event_table.index[event_table['Ligand Placed'] == True])[0])
                )
        except: 
            n_datasets_w_hit = 0

        try:    
            n_sites_w_hit = len(
                set(event_table['site_num'][event_table['Ligand Placed'] == True])
                )
        except: 
            n_sites_w_hit = 0

        # Hits per site

        n_hits_per_site = {
            i_site : sum(
                event_table["Ligand Placed"][event_table['site_num']==i_site]
                )
            for i_site in range(1, n_sites+1)
        }

        ###

        s_dict = dict(
            n_blobs = n_blobs,
            n_sites = n_sites,
            n_datasets = n_datasets,
            n_fitted = n_fitted,
            n_viewed = n_viewed,
            n_empty = n_empty,
            n_unviewed = n_unviewed,
            n_interesting_unfitted = n_interesting_unfitted,
            n_high_confidence = n_high_confidence,
            n_medium_confidence = n_medium_confidence,
            n_low_confidence = n_low_confidence,
            n_unranked = n_unranked,
            n_datasets_w_hit = n_datasets_w_hit,
            n_sites_w_hit = n_sites_w_hit,
            n_hits_per_site = n_hits_per_site,
        )

        return s_dict

    def make_table_block(self, title, table, width=12):

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

        header_title = 'PanDDA Inspection Summary'

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

