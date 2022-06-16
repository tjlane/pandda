import giant.logs as lg
logger = lg.getLogger(__name__)

import re
import pathlib as pl
import pandas as pd

from giant.mulch.labelling import (
    PathLabeller,
    )


class GetMetaFromTable(object):

    def __init__(self, meta_path):

        self.meta_path = pl.Path(meta_path)

        self.meta_table = pd.read_csv(
            str(self.meta_path)
            )

    def __call__(self, label):

        meta = self.meta_table.loc[label].todict()

        return meta


class GetMetaFromRegex(object):

    def __init__(self, meta_regex):

        self.regex_string = str(
            meta_regex
            )

        self.regex_compile = re.compile(
            str(meta_regex)
            )

    def __call__(self, string):

        m = self.regex_compile.search(
            str(string)
            )

        if m is None: 

            raise Exception(
                "Could not extract meta from \n\t{string}\n\t using regex '{regex}'".format(
                    string = string, 
                    regex = self.regex_string,
                    )
                )

        meta = m.groupdict()
        
        return meta


class ImportDatasets(object):

    def __init__(self,
        file_system,
        get_label = None,
        get_meta = None,
        use_meta_as_label = True, # TODO
        ):
        
        self.file_system = (
            file_system
            )

        self.get_label = (
            get_label
            )

        self.get_meta = (
            get_meta 
            if get_meta is not None
            else None
            )

        self.use_meta_as_label = bool(
            use_meta_as_label
            )

    def run(self, 
        data_paths,
        ):

        self.file_system.update()

        d_fs = self.file_system.datasets

        import_tuples = []

        logger.subheading('Validating datasets to be added')

        for m in data_paths: 

            logger(
                "MTZ: {mtz}".format(
                    mtz = m,
                    )
                )

            label = self.get_label(m)

            logger(
                "Label: {l}".format(
                    l = label,
                    )
                )

            if self.get_meta is not None: 

                meta_dict = self.get_meta(
                    m,
                    )

                logger(
                    'Meta: {d_string}\n'.format(
                        d_string = ', '.join([
                            "{k}={v}".format(k=k,v=v) 
                            for k,v in sorted(meta_dict.items())
                            ])
                        )
                    )

            else:

                meta_dict = None

            if self.use_meta_as_label is True: 
                
                assert meta_dict is not None
                
                label = '_'.join([
                    "{k}_{v}".format(k=k,v=v)
                    for k,v in sorted(
                        meta_dict.items()
                        )
                    ])

            assert not d_fs.exists(label), 'dataset already exists: {}'.format(
                label,
                )

            import_tuples.append(
                (label, m, meta_dict)
                )
    
        logger.subheading('Importing datasets')

        # now do the imports
        for label, m, meta_dict in import_tuples: 

            # add label to meta dict for simplicity
            if meta_dict is None: 
                meta_dict = {}

            if ('label' not in meta_dict): 
                meta_dict['label'] = label

            logger(
                "Importing {mtz} as {label}".format(
                    mtz = m,
                    label = label, 
                    )
                )

            dataset_folder = d_fs.create_dataset(
                label = label, 
                data_path = m,
                meta_dict = meta_dict,
                )

            logger(str(dataset_folder))

        # Check integrity of filesystem
        # d_fs.validate()

        return self.file_system
