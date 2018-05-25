import os, datetime


class FileManager(object):


    def __init__(self, rootdir):
        """Creates and stores output files to be neatly created and retrieved"""

        self._dir_parents = {'root':None}
        self._dir_names = {'root':rootdir}
        self._file_parents = {}
        self._file_names = {}

    def add_dir(self, dir_name, dir_tag, top_dir_tag=None, create=True, exists=True):
        """Store a directory `dir_name` under `dir_tag` in directory under `top_dir_tag`"""
        # Check that it hasn't been created already
        assert dir_tag not in self._dir_names.keys(), 'Directory already added: {}'.format(dir_tag)
        assert dir_tag not in self._dir_parents.keys(), 'Directory already added: {}'.format(dir_tag)
        # If no directory given, put in root
        if top_dir_tag is None: top_dir_tag = 'root'
        # Record the directory's parent and name
        self._dir_parents[dir_tag] = top_dir_tag
        self._dir_names[dir_tag] = dir_name
        # Create if it doesn't exist
        if create:
            self._make_directory_if_necessary(dir_tag)
        dir_path = self.get_dir(dir_tag)
        # Check it exists
        if exists:
            assert os.path.exists(dir_path)
        return dir_path

    def add_file(self, file_name, file_tag, dir_tag=None):
        """Store a filename `file_name` under `file_tag` in directory under `dir_tag`"""
        # Check that it hasn't been created already
        assert file_tag not in self._file_names.keys(), 'File already added: {}'.format(file_tag)
        assert file_tag not in self._file_parents.keys(), 'File already added: {}'.format(file_tag)
        # If no directory given, put in root
        if dir_tag is None: dir_tag = 'root'
        assert dir_tag in self._dir_names.keys(), 'Target directory does not exist: {}'.format(dir_tag)
        # Record the directory's parent and name
        self._file_parents[file_tag] = dir_tag
        self._file_names[file_tag] = file_name
        # Return filename
        return self.get_file(file_tag=file_tag)

    def get_dir(self, dir_tag):
        """Retrieve a dirname by it's dir_tag"""
        assert dir_tag in self._dir_names.keys(), 'Directory has not been added: {}'.format(dir_tag)
        assert dir_tag in self._dir_parents.keys(), 'Directory has not been added: {}'.format(dir_tag)
        # Extract the parents of the directory
        current = dir_tag
        path_segs = []
        while current is not None:
            path_segs.insert(0, self._dir_names[current])
            current = self._dir_parents[current]
        return os.path.join(*path_segs)

    def get_file(self, file_tag):
        """Retrieve a filename by it's file_tag"""
        assert file_tag in self._file_names.keys(), 'File has not been added: {}'.format(file_tag)
        assert file_tag in self._file_parents.keys(), 'File has not been added: {}'.format(file_tag)
        dir_tag = self._file_parents[file_tag]
        dir_path = self.get_dir(dir_tag=dir_tag)
        file_name = self._file_names[file_tag]
        return os.path.join(dir_path, file_name)

    def get_file_object(self, file_tag):
        """Retrieve a file object by it's file_tag"""
        assert file_tag in self._file_names.keys(), 'File has not been added: {}'.format(file_tag)
        assert file_tag in self._file_parents.keys(), 'File has not been added: {}'.format(file_tag)
        return FileObj(file=self.get_file(file_tag=file_tag), tag=file_tag)

    def check_and_create_directories(self):
        """Check that all directories exist and create them where necessary"""
        for dir_tag in self._dir_names.keys():
            self._make_directory_if_necessary(dir_tag=dir_tag)

    def _make_directory_if_necessary(self, dir_tag):
        dir_path = self.get_dir(dir_tag=dir_tag)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

    def change_root(self, new_root):
        self = copy.copy(self)
        self._dir_names['root'] = new_root
        return self


class FileObj(object):


    def __init__(self, file, tag=None):
        """File information object"""

        self.input = file

        self.path = os.path.abspath(file)
        self.tag = tag
        self.dir, self.name = os.path.split(self.path)
        self.base, self.ext = os.path.splitext(self.name)

    def __str__(self):
        return self.path

    def __call__(self):
        return self.path

    def get_datestamp(self, date='created'):
        stats = os.lstat(self.path)
        if date=='created':  return datetime.datetime.fromtimestamp(stats.st_ctime)
        if date=='modified': return datetime.datetime.fromtimestamp(stats.st_mtime)
        if date=='accessed': return datetime.datetime.fromtimestamp(stats.st_atime)
        raise Exception("Invalid option: date='{}'".format(date))


def compress_file(filename, delete_original=True):
    """Compress a file with gzip"""
    import gzip
    zip_file = filename + '.gz'
    with open(filename, 'r') as fh:
        f = gzip.open(zip_file, 'wb')
        f.write(fh.read())
        f.close()
    if delete_original:
        os.remove(filename)
    return zip_file

def decompress_all_files_in_directory(path, recursive=True, extension='.gz', delete_original=True):
    """Decompress in-place all files in a directory with a particular extension"""
    assert extension=='.gz', 'other extensions not yet implemented'
    assert os.path.exists(path), 'path does not exist: {}'.format(path)
    for cur_dir, _, cur_files in os.walk(path):
        for item in cur_files:
            if not item.endswith(extension): continue
            f_zip = os.path.join(cur_dir, item)
            f_out = os.path.splitext(f_zip)[0]
            assert not os.path.exists(f_out), 'output file already exists: {}'.format(f_out)
            if extension == '.gz':
                import gzip
                with gzip.open(f_zip) as fh_zip:
                    with open(f_out, 'w') as fh_out:
                        fh_out.write(fh_zip.read())
            else:
                pass
            assert os.path.exists(f_out), 'output file does not exist: {}'.format(f_out)
            if delete_original is True:
                os.remove(f_zip)
            #zip_obj = zipfile.ZipFile(item_path, 'r')
            #zip_obj.extractall(output_path)
            #zip_obj.close()
        # Breaking in top directory == non-recursive
        if recursive is False:
            break
