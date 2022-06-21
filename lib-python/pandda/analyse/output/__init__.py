import giant.logs as lg
logger = lg.getLogger(__name__)

import copy, collections
import pathlib as pl

from pandda.utils import (
    show_dict,
    merge_dicts,
    )

from .assign_sites import (
    AssignEventsToSites,
    )

from .statistics import (
    ExtractPanddaDatasetOutputInfo,
    )

from .tables import (
    MakePanddaResultsTables,
    )

from .graphs import (
    MakePanddaResultsGraphs,
    )

from .pymol import (
    MakePymolOutputImages,
    )

from .html import (
    MakePanddaResultsHtml,
    )


class PanddaResultsOutputter(object):
    """Output final results (outside processing loop)"""

    def __init__(self,
        output_dir,
        assign_sites = None,
        extract_info = None,
        write_tables = None,
        write_graphs = None,
        write_pymol_images = None,
        write_html = None,
        ):

        output_dir = pl.Path(output_dir)

        if (assign_sites is None):
            assign_sites = AssignEventsToSites()

        if (extract_info is None):
            extract_info = ExtractPanddaDatasetOutputInfo()

        if (write_tables is None):
            write_tables = MakePanddaResultsTables(
                output_directory = (
                    output_dir
                    ),
                )

        if (write_graphs is None):
            write_graphs = MakePanddaResultsGraphs(
                output_directory = (
                    output_dir / "graphs"
                    ),
                )

        if (write_pymol_images is None):
            write_pymol_images = MakePymolOutputImages(
                output_directory = (
                    output_dir
                    ),
                )

        if (write_html is None):
            write_html = MakePanddaResultsHtml(
                output_directory = (
                    output_dir
                    ),
                )

        self.output_dir = output_dir
        self.assign_sites = assign_sites
        self.extract_info = extract_info
        self.write_tables = write_tables
        self.write_graphs = write_graphs
        self.write_pymol_images = write_pymol_images
        self.write_html = write_html

    def __call__(self,
        mcd,
        pandda_results,
        ):

        if not self.output_dir.exists():
            self.output_dir.mkdir(parents=True)

        ###

        # Unpack input
        event_dicts = pandda_results['events']
        shell_dicts = pandda_results['shell_records']
        dataset_dicts = pandda_results['dataset_records']
        input_files = pandda_results['output_files']

        ###

        # Assign events to sites
        site_dicts = self.assign_sites(
            dataset_event_dicts = event_dicts, 
            )

        # Sort event dicts
        event_dicts = sorted(
            event_dicts,
            key = lambda e: (
                +e['site_idx'],
                -e['z_peak'],   # This is a key sorting point
                ),
            reverse = False,
            )

        ###

        output_files = collections.OrderedDict()

        logger.heading('Writing Pandda Output Files')

        dataset_dicts = self.extract_info(
            mcd = mcd,
            shell_dicts = shell_dicts,
            dataset_dicts = dataset_dicts, # gets overwritten with output variable
            )

        ###

        of = self.write_tables(
            dataset_dicts = dataset_dicts,
            event_dicts = event_dicts,
            site_dicts = site_dicts,
            )

        merge_dicts(
            master_dict = output_files,
            merge_dict = {'tables' : of},
            )

        ###

        of = self.write_graphs(
            event_dicts = event_dicts,
            shell_dicts = shell_dicts,
            dataset_dicts = dataset_dicts,
            )

        merge_dicts(
            master_dict = output_files,
            merge_dict = {'graphs' : of},
            )

        ###

        of = self.write_pymol_images(
            reference_structure = input_files['reference_files']['reference_model'],
            event_dicts = event_dicts,
            site_dicts = site_dicts,
            )

        merge_dicts(
            master_dict = output_files,
            merge_dict = {'graphs' : of},
            )

        ###

        combined_files = {}
        merge_dicts(master_dict=combined_files, merge_dict=input_files)
        merge_dicts(master_dict=combined_files, merge_dict=output_files)

        of = self.write_html(
            shell_dicts = shell_dicts,
            event_dicts = event_dicts,
            dataset_dicts = dataset_dicts,
            output_files = combined_files,
            )

        merge_dicts(
            master_dict = output_files,
            merge_dict = {'html' : of},
            )

        ###

        logger.subheading('Output Pandda Files')
        show_dict(output_files, logger=logger)

        return output_files
