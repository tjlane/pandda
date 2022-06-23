import giant.logs as lg
logger = lg.getLogger(__name__)

import json, collections
from libtbx import adopt_init_args
from libtbx.utils import Sorry, Failure


class JsonDataManager(object):


    def __init__(self,
        model_data,
        ):
        adopt_init_args(self, locals())
        self.validate()

    @classmethod
    def from_model_object(cls,
        model_object,
        ):
        """Extract json data from model object"""

        mo = model_object

        # Global dict
        model_data = collections.OrderedDict()

        # Level data
        level_data = collections.OrderedDict()
        model_data['level_descriptions'] = cls.extract_level_meta(mo)

        # TLS data
        tls_data = collections.OrderedDict()
        model_data['tls_level_data'] = cls.extract_tls_level_data(mo)

        return cls(
            model_data = model_data,
            )

    @classmethod
    def from_json_file(cls,
        filename,
        ):
        json_string = open(filename, 'r').read()
        return cls.from_json(
            json_string = json_string,
            )

    @classmethod
    def from_json(cls,
        json_string,
        ):
        model_data = json.loads(json_string)
        return cls(
            model_data = model_data,
            )

    def validate(self):

        for level_info in self.model_data.get('level_descriptions'):

            level_number = level_info.get('level_number')
            level_name = level_info.get('level_name')
            level_type = level_info.get('level_type')

            assert level_number is not None
            assert level_name is not None
            assert level_type is not None

            if level_type == 'tls':
                level_data = self.model_data.get('tls_level_data').get(level_name)
                assert level_data is not None
                assert level_data.get('level_number') == level_number
                assert level_data.get('level_name') == level_name
                self.validate_tls_level(level_data)
            elif level_type == 'adp':
                pass
            else:
                raise Exception('not implemented')

        # Check that all tls levels are included in the global list
        for tls_name in self.model_data.get('tls_level_data').keys():
            found = False
            for level_desc in self.model_data.get('level_descriptions'):
                if level_desc.get('level_type') != 'tls':
                    continue
                if level_desc.get('level_name') == tls_name:
                    found = True
                    break
            if found is False:
                raise Sorry('TLS Level "{}" is not found in the level_descriptions data'.format(tls_name))

        return

    def validate_tls_level(self, data):
        """Check the fields for a tls level"""

        required_fields = ['level_number', 'level_name', 'tls_group_data']
        required_types  = [int, str, list]

        for r_field, r_type in zip(required_fields, required_types):
            assert r_field in data
            assert isinstance(data.get(r_field), r_type)

        for group_data in data.get('tls_group_data'):
            self.validate_tls_group(group_data)

    def validate_tls_group(self, data):
        """Check the fields for a tls group"""
        required_fields = ['group_number', 'selection', 'tls_origins', 'tls_modes']
        required_types  = [int, str, dict, list]
        for r_field, r_type in zip(required_fields, required_types):
            assert r_field in data
            assert isinstance(data.get(r_field), r_type)

        for mode_data in data.get('tls_modes'):
            self.validate_tls_mode(mode_data)

    def validate_tls_mode(self, data):
        """Check the fields for a tls mode"""
        required_fields = ['T', 'L', 'S', 'amplitudes']
        required_types  = [list, list, list, dict]

        for r_field, r_type in zip(required_fields, required_types):

            assert r_field in data
            assert isinstance(data.get(r_field), r_type)

        assert len(data.get('T')) == 6
        assert len(data.get('L')) == 6
        assert len(data.get('S')) == 9

    @staticmethod
    def extract_level_meta(model_object):

        mo = model_object

        all_level_meta = []

        for i_l in range(len(mo.all_level_names)):

            level_meta = collections.OrderedDict()
            all_level_meta.append(level_meta)

            level_meta['level_number'] = i_l+1
            level_meta['level_name'] = str(model_object.all_level_names[i_l])
            level_meta['level_type'] = str(model_object.all_level_types[i_l])

        return all_level_meta

    @staticmethod
    def extract_tls_level_data(model_object):

        mo = model_object

        dataset_labels = list(map(str, mo.dataset_labels))

        tls_level_data = collections.OrderedDict()

        for i_l, level_name in enumerate(mo.all_level_names):

            # Only select tls levels
            if mo.all_level_types[i_l] != 'tls':
                continue

            # Extract the index of this relative to the tls levels
            i_tls_level = mo.tls_level_names.index(level_name)

            level_data = collections.OrderedDict()
            tls_level_data[str(level_name)] = level_data

            level_data['level_number']   = i_l+1
            level_data['level_name']     = str(level_name)
            level_data['tls_group_data'] = []

            tls_strings = mo.tls_selection_strings[i_tls_level]
            tls_objects = mo.tls_objects[i_tls_level]

            # Iterate through groups in level
            for i_group, (tls_s, tls_g) in enumerate(zip(tls_strings, tls_objects)):

                group_data = collections.OrderedDict()
                level_data['tls_group_data'].append(group_data)

                group_data['group_number']    = i_group+1
                group_data['selection']       = str(tls_s)
                group_data['number_of_atoms'] = tls_g.n_atoms
                group_data['tls_modes']       = []
                group_data['tls_origins']     = collections.OrderedDict(zip(dataset_labels, tls_g.origins))

                for i_mode, ma_values in enumerate(tls_g.tls_parameters):

                    mode_data = collections.OrderedDict()
                    group_data['tls_modes'].append(mode_data)

                    mode_data['index'] = i_mode

                    for letter in 'TLS':
                        matrix_values = list(ma_values.matrices.get(letter))
                        mode_data[str(letter)] = matrix_values

                    amplitudes = list(ma_values.amplitudes.get())
                    amplitudes_dict = collections.OrderedDict(zip(dataset_labels, amplitudes))
                    mode_data['amplitudes'] = amplitudes_dict

        return tls_level_data

    def as_json(self):
        return json.dumps(self.model_data, indent=2)

    def write_json(self, filename, mode='a'):
        json_string = self.as_json()
        with open(filename, mode=mode) as fh:
            fh.write(json_string)

    def apply_to_model_object(self, model_object):
        """Insert TLS values into a model object"""

        mo = model_object

        msg = """Non-critical warning while applying model data\n\t> {}"""

        # Datasets that are not in the new model (the one being updated)
        missing_datasets = set()

        # Apply the tls levels
        for level_name, level_data in self.model_data.get('tls_level_data').items():

            if level_name not in mo.tls_level_names:
                logger.warning(msg.format('Level "{}" not found in model'.format(level_name)))
                continue

            # Find the index of this level
            i_tls_level = mo.tls_level_names.index(level_name)

            # Extract selection strings and object
            tls_selection_strings = mo.tls_selection_strings[i_tls_level]
            tls_objects = mo.tls_objects[i_tls_level]

            # Iterate through and insert information
            for source_info in level_data.get('tls_group_data'):
                # Selection string from json
                sel_string = source_info.get('selection')

                if sel_string not in tls_selection_strings:
                    logger.warning(msg.format('Group "{}" not found in selection strings for level "{}"'.format(sel_string, level_name)))
                    continue

                # Find the index of the tls group
                tls_i = tls_selection_strings.index(sel_string)

                # Extract the group object
                tls_g = tls_objects[tls_i]

                # Apply matrices and amplitudes
                for mode_info in source_info.get('tls_modes'):
                    # Index of this mode in the new model
                    i_mode = mode_info.get('index')
                    # Error if longer
                    if i_mode >= tls_g.n_modes:
                        raise Sorry('Number of TLS modes ({}) in JSON data is greater than in the model.'.format(i_mode+1, tls_g.n_modes))

                    # tls mode object
                    tls_values = tls_g.tls_parameters[i_mode]

                    # Apply tls matrices values
                    tls_values.matrices.set(values=mode_info.get('T'), component_string='T')
                    tls_values.matrices.set(values=mode_info.get('L'), component_string='L')
                    tls_values.matrices.set(values=mode_info.get('S'), component_string='S')

                    # Apply tls amplitudes values
                    for dataset_label, amplitude in mode_info.get('amplitudes').items():
                        # Check that this dataset exists in the new model
                        if (dataset_label not in model_object.dataset_labels):
                            skip_datasets.add(dataset_label)
                            continue
                        # Extract index of this dataset in the new model
                        i_dst = mo.dataset_labels.index(dataset_label)
                        # Apply amplitudes
                        tls_values.amplitudes.set(values=[amplitude], selection=[i_dst])

                # Apply origins
                for dataset_label, origin in source_info.get('tls_origins').items():
                    # Check that this dataset exists in the new model
                    if (dataset_label not in model_object.dataset_labels):
                        skip_datasets.add(dataset_label)
                        continue
                    # Extract index of this dataset in the new model
                    i_dst = mo.dataset_labels.index(dataset_label)
                    # Apply origin
                    tls_g.origins[i_dst] = origin

        if len(missing_datasets) > 0:
            txt = "The following datasets are present in the JSON data but not in the new model: \n\t\t{}".format('\n\t\t'.join(missing_datasets))
            logger.warning(msg.format(txt))


def write_model_as_json(model_object, output_file):

    output_json_manager = JsonDataManager.from_model_object(
        model_object = model_object,
        )
    output_json_manager.write_json(
        filename = output_file,
        mode = 'w',
        )
