import copy, collections
from libtbx import adopt_init_args


class WriteNewCustomLevelsEffFile(object):


    def __init__(self,
        master_phil,
        custom_level_scope_name = 'model.custom_level',
        ):
        adopt_init_args(self, locals())

        # Extract the required scope
        self.custom_level_scope = self.master_phil.get(custom_level_scope_name).copy()
        self.custom_level_scope.name = '.'.join(custom_level_scope_name.split('.')[:-1]) # hack -- dangerous?

    def __call__(self,
        selection_string_hierarchy,
        output_filename,
        output_level_names,
        comment_lines = [],
        depth = None,
        insert_before = None,
        insert_after = None,
        reverse_levels = False,
        mode = 'a',
        ):
        """Take levels of selection strings and convert to phil strings"""

        assert mode in ['a', 'w']

        # Python object for custom_level scope (list of template item)
        custom_level_scope_extract = self.custom_level_scope.extract()
        # Python object for each custom level
        custom_level_object = custom_level_scope_extract.custom_level.pop()

        assert len(custom_level_scope_extract.custom_level) == 0

        custom_level_list = []
        for i_l, string_set in enumerate(selection_string_hierarchy):
            if not string_set: continue
            custom_level_new = copy.copy(custom_level_object)
            custom_level_new.depth = depth
            custom_level_new.insert_before = insert_before
            custom_level_new.insert_after = insert_after
            custom_level_new.label = output_level_names[i_l]
            custom_level_new.selection = string_set
            custom_level_list.append(custom_level_new)

        if reverse_levels is True:
            custom_level_list.reverse()

        # Put custom level list back in
        custom_level_scope_extract.custom_level = custom_level_list

        custom_level_formatted = self.custom_level_scope.format(custom_level_scope_extract)

        with open(output_filename, mode) as fh:
            for c in comment_lines:
                fh.write('# '+c+'\n')
            fh.write(custom_level_formatted.as_str())

        return custom_level_formatted


class _MakeNewEffFile(object):


    def __init__(self,
        selection_strings_dict,
        master_phil,
        custom_level_scope_name = 'model.custom_level',
        ):

        self.make_eff_file = WriteNewCustomLevelsEffFile(
            master_phil = master_phil,
            custom_level_scope_name = custom_level_scope_name,
            )

        adopt_init_args(self, locals())

    @classmethod
    def from_model_object(cls,
        model_object,
        master_phil,
        custom_level_scope_name = 'model.custom_level',
        ):

        # to be removed once tls_selection_strings is stored as a dict
        assert len(model_object.tls_level_names) == len(model_object.tls_selection_strings)
        selection_strings_dict = dict(zip(model_object.tls_level_names, model_object.tls_selection_strings))

        return cls(
            selection_strings_dict = selection_strings_dict,
            master_phil = master_phil,
            custom_level_scope_name = custom_level_scope_name,
            )


class MakeNewCustomLevelEffFilesFromIndices(_MakeNewEffFile):


    def __call__(self,
        level_name,
        indices_hierarchy,
        output_filename,
        label_template,
        comment_lines = [],
        depth = None,
        insert_before = None,
        insert_after = None,
        reverse_levels = False,
        ):

        # Skip missing levels (allowed)
        if level_name not in self.selection_strings_dict:
            valid_keys = ', '.join(map(str,self.selection_strings_dict.keys()))
            raise Exception('level_name "{}" not found in selection_strings_dict (valid keys: {})'.format(level_name, valid_keys))

        # Collections of strings for each level
        strings_hierarchy = self.convert_indices_hierarchy_to_selection_strings(
            indices_hierarchy = indices_hierarchy,
            selection_strings = self.selection_strings_dict[level_name],
            )

        # level1: '1+2,4+5+6,...'
        output_level_name_strings = [', '.join(['+'.join(map(str,[i+1 for i in g_indices])) for g_indices in l_indices]) for l_indices in indices_hierarchy]

        eff_phil = self.make_eff_file(
            selection_string_hierarchy = strings_hierarchy,
            output_filename = output_filename,
            output_level_names = list(map(label_template.format, output_level_name_strings)),#[label_template.format(i_l+1) for i_l in range(len(strings_hierarchy))],
            comment_lines = comment_lines,
            depth = depth,
            insert_before = insert_before,
            insert_after = insert_after,
            reverse_levels = reverse_levels,
            )

        return eff_phil.as_str()

    @staticmethod
    def convert_indices_hierarchy_to_selection_strings(
        indices_hierarchy,
        selection_strings,
        ):
        h_strings = []
        for level_sets in indices_hierarchy:
            l_strings = []
            for group in level_sets:
                g_string = '('+') or ('.join([selection_strings[i] for i in group])+')'
                l_strings.append(g_string)
            h_strings.append(tuple(l_strings))
        return h_strings


class MakeNewCustomHierarchyEffFilesFromIndices(_MakeNewEffFile):


    def __call__(self,
        indices_hierarchy,
        input_level_names,
        output_level_names,
        output_filename,
        comment_lines = [],
        depth = None,
        insert_before = None,
        insert_after = None,
        reverse_levels = False,
        ):

        # Collections of strings for each level
        strings_hierarchy = self.convert_indices_hierarchy_to_selection_strings(
            indices_hierarchy = indices_hierarchy,
            selection_strings = [self.selection_strings_dict[level_name] for level_name in input_level_names],
            )

        eff_phil = self.make_eff_file(
            selection_string_hierarchy = strings_hierarchy,
            output_filename = output_filename,
            output_level_names = output_level_names,
            comment_lines = comment_lines,
            depth = depth,
            insert_before = insert_before,
            insert_after = insert_after,
            reverse_levels = reverse_levels,
            )

        return eff_phil.as_str()

    @staticmethod
    def convert_indices_hierarchy_to_selection_strings(
        indices_hierarchy,
        selection_strings,
        ):
        h_strings = []
        for level_sets, level_strings in zip(indices_hierarchy, selection_strings):
            l_strings = []
            for group in level_sets:
                g_string = '('+') or ('.join([level_strings[i] for i in group])+')'
                l_strings.append(g_string)
            h_strings.append(tuple(l_strings))
        return h_strings

