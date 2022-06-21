import json 

import pathlib as pl


class DumpDictToJson(object):

    def __init__(self, output_path):

        self.output_path = pl.Path(output_path)

    def __call__(self, records_dict):

        self.write_json(records=records_dict)

        return self.output_path

    def write_json(self, records):

        json_string = json.dumps(records, indent=2)

        with open(str(self.output_path), "w") as f:
            f.write(json_string)
