import json
import multiprocessing
import pathlib as pl

from giant.processors import (
    Processor,
    ProcessorJoblib,
    )

system_cpus = multiprocessing.cpu_count()


class UpdateConfiguration(object):

    def __init__(self,
        backend,
        ):

        self.backend = backend

    def __call__(self, configuration):

        configuration.data['backend'] = self.backend

        configuration.refresh() # check it works

        configuration.write() # then write

        return configuration


class MassRefineConfiguration(object): 

    def __init__(self, json_path):

        self.json_path = pl.Path(json_path)

        self.data = self.read()

        self.refresh()

    def read(self):

        if not self.json_path.exists():
            
            return {}

        with open(str(self.json_path), 'r') as fh:

            data = json.loads(
                fh.read()
                )

        return data

    def write(self):

        with open(str(self.json_path), 'w') as fh:

            fh.write(
                json.dumps(self.data)
                )

    def refresh(self):

        self.processor = self.get_processor()

    def get_processor(self):

        d = self.data

        if d.get('backend') == "parallel_joblib":
            proc = ProcessorJoblib(
                n_cpus = d.get("n_cpus", system_cpus),
                )
        elif d.get('backend') == "serial": 
            proc = Processor()
        else:
            proc = Processor()

        return proc

