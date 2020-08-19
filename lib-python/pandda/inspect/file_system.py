import pathlib as pl


class PanddaInspectorFilenames:

    def __init__(self, 
        input_directory,
        output_directory,
        ):

        input_directory = pl.Path(
            input_directory
            )

        self.input_event_csv = str(
            input_directory / "results" / "pandda_analyse_events.csv"
            )

        self.input_site_csv = str(
            input_directory / "results" / "pandda_analyse_sites.csv"
            )

        ###

        output_directory = pl.Path(
            output_directory
            )

        self.output_event_csv = str(
            output_directory / "inspect" / "pandda_inspect_events.csv"
            )

        self.output_site_csv = str(
            output_directory / "inspect" / "pandda_inspect_sites.csv"
            )