import json
import pathlib as pl


class FolderWithMeta(object): 

    meta_path = None

    def get_meta(self):

        json_string = open(str(self.meta_path),'r').read()

        return json.loads(
            json_string
            )

    def write_meta(self, meta_dict):

        with open(str(self.meta_path), 'w') as fh: 
            fh.write(
                json.dumps(
                    meta_dict
                    )
                )

    def update_meta(self, update_dict):
        meta = self.get_meta()
        meta.update(update_dict)
        self.write_meta(meta)
        return meta

    def tmp_dir(self):

        tmp_dir = (
            self.root_dir / "tmp_files"
            )

        if not tmp_dir.exists():
            tmp_dir.mkdir()

        return tmp_dir


class DirectoryPathGenerator(object):

    int_format = "{:04d}"
    prefix = ""

    def __init__(self, root_dir):

        self.root_dir = pl.Path(root_dir)

        assert '_' not in self.int_format, self.int_format
        assert '_' not in self.prefix, self.prefix

    def get_next(self, suffix=None):

        if suffix is not None:

            suffix = str(suffix).replace('/','_')

            # if not suffix.startswith('_'): 
            #     suffix = ('_'+suffix)

        else: 

            suffix = ''

        next_int = self.get_next_int()

        new_path = self.root_dir / (
            str(self.prefix) + '_' + 
            str(self.int_format).format(next_int) + '_' +
            str(suffix)
            )

        return new_path

    def get_all(self):

        return [
            p for p in sorted(
                self.root_dir.glob(
                    self.prefix+'*',
                    )
                )
            if p.is_dir()
            ]

    def get_all_as_dict_by_int(self):

        return {
            self.get_int_from_name(p.name) : p
            for p in self.get_all()
            }
    
    def get_all_as_dict_by_name(self):

        return {
            self.get_label_from_name(p.name) : p
            for p in self.get_all()
            }

    def get_next_int(self):

        return self.get_max_int() + 1

    def get_max_int(self): 

        all_ints = []

        for p in self.get_all():

            try: 
                i = self.get_int_from_name(p.name)
            except: 
                continue

            all_ints.append(i)

        if len(all_ints) == 0: 
            return 0
            
        return max(all_ints)

    def get_int_from_name(self, name):

        return int(
            name.split('_')[1] # 0prefix_1int_2suffix
            )

    def get_label_from_name(self, name):

        return '_'.join(
            name.split('_')[2:] # 0prefix_1int_2suffix
            )
