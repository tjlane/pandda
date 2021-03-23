import json
import pathlib as pl


class GetPanddaInspectInputOutputFiles:

    def __init__(self, 
        input_directory,
        output_directory,
        ):

        self.input_directory = pl.Path(
            input_directory
            )

        self.output_directory = pl.Path(
            output_directory
            )

        assert self.input_directory.exists()

        if not self.output_directory.exists():
            self.output_directory.mkdir(parents=True)

    def __call__(self):

        results_dict = self.get_json_data()

        output_files = results_dict['output_files']

        ret_dict = {
            'pandda_files' : output_files,
            }

        ret_dict.update(
            self.find_input_files(output_files)
            )

        ret_dict.update(
            self.setup_output_files()
            )

        return ret_dict

    def get_json_data(self):

        results_json = (
            self.input_directory / "results.json"
            )

        results_dict = json.loads(
            open(str(results_json), 'r').read()
            )

        return results_dict

    def find_input_files(self, output_files):

        input_event_csv = output_files['tables']['events_table']

        remove_prefix = self.find_output_prefix(input_event_csv)

        input_events_csv = (
            self.input_directory / (
                pl.Path(
                    input_event_csv
                    ).relative_to(remove_prefix)
                )
            )

        input_sites_csv = (
            self.input_directory / (
                pl.Path(
                    output_files['tables']['sites_table']
                    ).relative_to(remove_prefix)
                )
            )
        
        assert input_events_csv.exists()
        assert input_sites_csv.exists()

        return {
            'input_events' : input_events_csv,
            'input_sites' : input_sites_csv,
            'remove_prefix' : remove_prefix,
        }

    def setup_output_files(self):

        out_dir = (
            self.output_directory
            )

        if not out_dir.exists():
            out_dir.mkdir(parents=True)

        output_events_csv = str(
            out_dir / "pandda_inspect_events.csv"
            )

        output_sites_csv = str(
            out_dir / "pandda_inspect_sites.csv"
            )

        return {
            'output_events' : output_events_csv,
            'output_sites' : output_sites_csv,
        }

    def find_output_prefix(self, filepath):

        filepath = pl.Path(filepath)

        top_dir = self.input_directory

        prefix = filepath.parent

        for _ in filepath.parts:

            test_f = filepath.relative_to(prefix)

            if (top_dir / test_f).exists():
                break

            prefix = prefix.parent
        
        return prefix


