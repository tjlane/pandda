import giant.logs as lg
logger = lg.getLogger(__name__)

import os
import pathlib as pl

from giant.mulch.collection import (
    MultiCrystalDataset,
    CrystallographicDataset,
    )
from giant.mulch.checks import (
    CrystallographicDatasetChecker,
    )


class DefaultInputDirectoryParser(object):

    def __init__(self,
        data_dirs,
        pdb_style,
        mtz_style,
        pdb_regex = None,
        mtz_regex = None,
        dir_regex = None,
        only_datasets = None,
        ignore_datasets = None,
        dataset_prefix = None,
        ):

        self.name = "DefaultParseInputDirectory"

        self.data_dirs = data_dirs.rstrip('/') # just in case it has a leading /
        self.pdb_style = pdb_style.strip('/')
        self.mtz_style = mtz_style.strip('/')
        self.pdb_regex = pdb_regex
        self.mtz_regex = mtz_regex
        self.dir_regex = dir_regex
        self.only_datasets = only_datasets
        self.ignore_datasets = ignore_datasets
        self.dataset_prefix = dataset_prefix

        self.new_files = None
        self.filtered_files = None
        self.empty_directories = None

    def __call__(self):
        """Builds a list of input files from the command line arguments passed"""

        import glob

        # ==============================>
        # Find datasets in the input directories
        # ==============================>
        new_files = {}
        filtered_files = {}
        empty_directories = []

        for d_dir in sorted(glob.glob(self.data_dirs)):

            d_dir = pl.Path(d_dir)

            pdb_files = [p for p in d_dir.glob(self.pdb_style) if p.exists()]
            mtz_files = [p for p in d_dir.glob(self.mtz_style) if p.exists()]

            if (not pdb_files) or (not mtz_files):

                if (not pdb_files):
                    logger("No PDB file found in {d}".format(d=str(d_dir)))

                if (not mtz_files):
                    logger("No MTZ file found in {d}".format(d=str(d_dir)))

                empty_directories.append(str(d_dir))

                continue

            elif (len(pdb_files) > 1) or (len(mtz_files) > 1):

                if len(pdb_files) > 1:
                    logger(
                        'More than one matching PDB file found for {glob}: \n{files}'.format(
                            glob = str(d_dir / self.pdb_style),
                            files = str(pdb_files),
                            )
                        )

                if len(mtz_files) > 1: 
                    logger(
                        'More than one matching MTZ file found for {glob}: \n{files}'.format(
                            glob = str(d_dir / self.mtz_style),
                            files = str(pdb_files),
                            )
                        )

                raise ValueError(
                    "Multiple matching files found in input directory ({pdbs} and {mtzs})".format(
                        pdbs = str(pdb_files),
                        mtzs = str(mtz_files),
                        )
                    )

            # Found PDB and MTZ file in directory
            pdb_path = str(pdb_files[0])
            mtz_path = str(mtz_files[0])
    
            # Find tag from file names
            dataset_tag = self.find_dataset_tag(
                pdb_path = pdb_path,
                mtz_path = mtz_path,
                )

            # Check it doesn't already exist
            if (dataset_tag in new_files): 
                raise ValueError(
                    "Dataset with this label has already been found -- all datasets must have unique labels. \n" + \
                    "This dataset: label {}, pdb {}, mtz {}\n".format(
                        dataset_tag, pdb_path, mtz_path,
                        ) + \
                    "Previous dataset: label {}, pdb {}, mtz {}\n".format(
                        dataset_tag, new_files[dataset_tag][0], new_files[dataset_tag][1],
                        )
                    )

            new_files[dataset_tag] = (str(pdb_path), str(mtz_path))

        # Select only those requested
        if (self.only_datasets is not None):

            remove_tags = set(
                new_files.keys()
                ).difference(
                list(self.only_datasets)
                )

            for t in remove_tags: 
                filtered_files[t] = ("this dataset label is not in 'only_datasets'", new_files.pop(t))

            if len(remove_tags) == 0: 
                logger("No datasets were removed using the only_datasets filter")

            if len(new_files) == 0:
                raise ValueError(
                    "After filtering with only_datasets, no files are remaining.\n" + \
                    "Selected dataset labels: {}\n".format(list(self.only_datasets)) + \
                    "Removed dataset labels: {}\n".format(list(remove_tags))
                    )

        # Filter out manually labelled datasets to ignore
        if (self.ignore_datasets is not None):

            remove_tags = set(
                new_files.keys()
                ).intersection(
                list(self.ignore_datasets)
                )

            for t in remove_tags: 
                filtered_files[t] = ("this dataset label is in 'ignore_datasets'", new_files.pop(t))

            if len(remove_tags) == 0:
                logger("No datasets were removed using the ignore_datasets filter")

            if len(new_files) == 0: 
                raise ValueError(
                    "After filtering with ignore_datasets, no files are remaining.\n" + \
                    "Selected dataset labels: {}\n".format(list(self.ignore_datasets)) + \
                    "Removed dataset labels: {}\n".format(list(remove_tags))
                    )

        self.new_files = new_files
        self.filtered_files = filtered_files
        self.empty_directories = empty_directories

        return self.new_files

    def __str__(self):

        s_ = (
            'ParserType: {name}\n'
            '----------\n'
            'Datasets Found: {n_found}\n'
            'Datasets Filtered: {n_filtered}\n'
            'Empty Directories: {n_empty}\n'
            '----------\n'
            'Found Files: \n\t{found_files}\n'
            '----------\n'
            'Filtered Files: \n\t{filtered_files}\n'
            '----------\n'
            'Empty Directories: \n\t{empty_directories}\n'
            '----------\n'
            ).format(
            name = self.name,
            n_found = len(self.new_files),
            found_files = (
                '\n'.join([
                    '{d}:\n\t{p}\n\t{m}'.format(
                        d=dtag, p=pdb, m=mtz,
                        )
                    for dtag, (pdb, mtz) 
                    in sorted(self.new_files.items())
                    ])
                if len(self.new_files)
                else str(None)
                ).strip().replace('\n','\n\t'),
            n_filtered = len(self.filtered_files),
            filtered_files = (
                '\n'.join([
                    '{d}:\n\treason: {r},\n\t{p},\n\t{m}'.format(
                        d=dtag, r=reason, p=pdb, m=mtz,
                        )
                    for dtag, (reason, (pdb, mtz)) 
                    in sorted(self.filtered_files.items())
                    ])
                if len(self.filtered_files)
                else str(None)
                ).strip().replace('\n','\n\t'),
            n_empty = len(self.empty_directories),
            empty_directories = (
                '\n'.join(
                    sorted(self.empty_directories)
                    )
                if len(self.empty_directories)
                else str(None)
                ).strip().replace('\n','\n\t'),
            )

        return s_

    def info(self):

        return {
            'empty_directories' : self.empty_directories,
            'filtered_files' : self.filtered_files,
            'new_files' : self.new_files,
        }

    def find_dataset_tag(self, pdb_path, mtz_path):

        pdb_tag = mtz_tag = dir_tag = None

        if ('*' in self.pdb_style):

            pdb_tag = self.extract_regex(
                string = os.path.basename(pdb_path),
                regex = (
                    self.pdb_regex 
                    if (self.pdb_regex is not None)
                    else self.pdb_style.replace('*', '(.*)')
                    ),
                )

        if ('*' in self.mtz_style):

            mtz_tag = self.extract_regex(
                string = os.path.basename(mtz_path),
                regex = (
                    self.mtz_regex 
                    if (self.mtz_regex is not None)
                    else self.mtz_style.replace('*', '(.*)')
                    ),
                )
                
        if ('*' in self.data_dirs): 

            dir_tag = self.extract_regex(
                string = os.path.dirname(pdb_path),
                regex = (
                    self.dir_regex 
                    if (self.dir_regex is not None)
                    else self.data_dirs.replace('*', '(.*)')
                    ),
                )

        if [pdb_tag, mtz_tag, dir_tag].count(None) == 3: 
            raise ValueError("No tags were identified!")

        tag = self.select_dataset_tag(
            dir_tag = dir_tag,
            pdb_tag = pdb_tag,
            mtz_tag = mtz_tag,
            )

        assert (tag is not None)

        if (self.dataset_prefix is not None):
            tag = (self.dataset_prefix + tag)

        return tag

    def extract_regex(self, string, regex):

        import re

        matches = re.findall(regex, string)                    

        if len(matches) == 0: 
            raise ValueError(
                "No matches for regex {regex} in {string}".format(
                    regex = regex,
                    string = string,
                    )
                )

        tag_set = set()
        for m in matches:
            if isinstance(m, tuple):
                tag_set.update(m)
            else: 
                tag_set.add(m)

        # Multiple nonredundant tags -- error
        if len(tag_set) > 1: 
            raise ValueError(
                "Multiple matches for {regex} in {string}: {matches}".format(
                    regex = regex, 
                    string = string,
                    matches = list(tag_set),
                    )
                )

        # Only one tag -- extract from list
        tag = str(list(tag_set)[0])

        # Check for null string
        if len(tag) == 0: 
            tag = None

        return tag

    def select_dataset_tag(self, 
        dir_tag, 
        pdb_tag,
        mtz_tag,
        ):

        #####
        # Consistency checks! 
        #####
        differences = False
        if (dir_tag is not None) and (pdb_tag is not None) and (dir_tag != pdb_tag): 
            logger.warning(
                "Extracted labels from directory and pdb file are not the same: '{d}' != '{f}'".format(
                    d = dir_tag, f = pdb_tag,
                    )
                )
            differences = True

        if (dir_tag is not None) and (mtz_tag is not None) and (dir_tag != mtz_tag):
            logger.warning(
                "Extracted labels from directory and mtz file are not the same: '{d}' != '{f}'".format(
                    d = dir_tag, f = mtz_tag,
                    )
                )
            differences = True

        if (pdb_tag is not None) and (mtz_tag is not None) and (pdb_tag != mtz_tag):
            logger.warning(
                "Extracted labels from pdb file and mtz file are not the same: '{f1}' != '{f2}'".format(
                    f1 = pdb_tag, f2 = mtz_tag,
                    )
                )
            differences = True

        #####
        # Select a tag in a preferred order (regex first, then dir -> pdb -> mtz)
        #####
        if (dir_tag is not None) and (self.dir_regex is not None):

            tag = dir_tag
            tag_used = "directory (from dir_regex)"

        elif (pdb_tag is not None) and (self.pdb_regex is not None):

            tag = pdb_tag
            tag_used = "pdb file (from pdb_regex)"

        elif (mtz_tag is not None) and (self.mtz_regex is not None):

            tag = mtz_tag
            tag_used = "mtz file (from mtz_regex)"

        elif (dir_tag is not None):

            tag = dir_tag
            tag_used = "directory (from * in data_dirs)"

        elif (pdb_tag is not None):

            tag = pdb_tag
            tag_used = "pdb file (from * in pdb_style)"

        elif (mtz_tag is not None):

            tag = mtz_tag
            tag_used = "mtz file (from * in mtz_style)"

        if (differences is True):
            logger.warning(
                "Different tags were identified -- using {used} labelling: {tag}".format(
                    used = tag_used, tag = tag,
                    )
                )

        return tag


class DefaultDatasetLoader(object):

    DatasetClass = CrystallographicDataset
    DatasetCheckerClass = CrystallographicDatasetChecker

    def __init__(self):

        self.name = "DefaultDatasetLoader"

    def __call__(self, new_files):
        """Read in maps for the input datasets"""

        # ==============================>
        # Load datasets in parallel
        # ==============================>
        loaded_datasets = {
            dtag: self.DatasetClass.from_file(
                model_filename = pdb,
                data_filename = mtz,
            ).label(
                num = num,
                tag = dtag,
            )
            for num, (dtag, (pdb, mtz)) in enumerate(new_files.items())
        }

        checker = self.DatasetCheckerClass()

        errors = []
        for dtag, d in loaded_datasets.items():
            messages = checker(d)
            for m in messages:
                errors.append("{}: {}".format(dtag, m))

        if len(errors) > 0:
            for e in errors:
                logger.error(e)
            raise ValueError("Loaded datasets do not pass checks. See above.")

        return loaded_datasets

    def __str__(self):

        return ""


class MultiDatasetDataloader(object):

    InputDirectoryParser = DefaultInputDirectoryParser
    DatasetLoader = DefaultDatasetLoader

    def __init__(self,
        data_dirs,
        pdb_style,
        mtz_style,
        pdb_regex,
        mtz_regex,
        dir_regex,
        only_datasets,
        ignore_datasets,
        dataset_prefix,
        ):

        self.name = "DefaultDataloader"

        assert [pdb_style, mtz_style].count(None) < 2, 'Must provide at least one'

        if (mtz_style is None): 
            mtz_style = (os.path.splitext(pdb_style)[0] + '.mtz')

        if (pdb_style is None): 
            pdb_style = (os.path.splitext(mtz_style)[0] + '.pdb')

        self.parse_input_directories = self.InputDirectoryParser(
            data_dirs = data_dirs,
            pdb_style = pdb_style,
            mtz_style = mtz_style,
            pdb_regex = pdb_regex,
            mtz_regex = mtz_regex,
            dir_regex = dir_regex,
            only_datasets = only_datasets,
            ignore_datasets = ignore_datasets,
            dataset_prefix = dataset_prefix,
        )

        self.load_datasets = self.DatasetLoader()

    def __str__(self):

        s_ = (
            'DataLoaderType: {name}\n'
            'Parsed Directories:\n\t{parser}\n'
            'Loaded Datasets:\n\t{loader}\n'
            ).format(
            name = self.name,
            parser = str(
                self.parse_input_directories
                ).strip('\n').replace('\n','\n\t'), 
            loader = str(
                self.load_datasets
                ).strip('\n').replace('\n','\n\t'),
            )

        return s_.strip('\n')

    def __call__(self):

        file_triplets = self.parse_input_directories()

        datasets = self.load_datasets(file_triplets)

        mcd = MultiCrystalDataset(datasets)

        return mcd

    def info(self):

        info = {}
        info.update(self.parse_input_directories.info())

        return info
