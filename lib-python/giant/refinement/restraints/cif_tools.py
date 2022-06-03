import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy as np
import gemmi as gm
import re


class CifManager:

    LINK_LIST_NAME = 'link_list'
    MOD_LIST_NAME = 'mod_list'

    LINK_BLOCK_PRE = 'link_'
    MOD_BLOCK_PRE = 'mod_'

    CHEM_MOD_ROOT = '_chem_mod'
    CHEM_LINK_ROOT = '_chem_link'

    def __init__(self, cif_obj):

        self.cif_obj = cif_obj

    def __str__(self):

        block_string = '\n'.join([
            b.name for b in self
            ])

        list_block_string = '\n'.join([
            'Block {block_name}\n\t{block_string}'.format(
                block_name = b.name,
                block_string = (
                    b.as_string()
                    ).strip('\n').replace('\n','\n\t'),
                )
            for b in self
            if b.name.endswith('_list')
            ])

        s_ = (
            'Cif object\n'
            '| Source: {src}\n'
            '| Block List: \n'
            '|\t{block_string}\n'
            '| List Blocks: \n'
            '|\t{list_block_string}\n'
            '`---->'
            ).format(
            src = self.cif_obj.source,
            block_string = (
                block_string
                ).strip('\n').replace('\n','\n|\t'),
            list_block_string = (
                list_block_string
                ).strip('\n').replace('\n','\n|\t'),
            )

        return s_.strip()

    def __iter__(self):

        for a in self.cif_obj:
            yield a

    @classmethod
    def from_file(cls, filepath):

        cif_obj = gm.cif.read(
            str(filepath)
            )

        return cls(cif_obj=cif_obj)

    def write_cif(self, filepath):

        self.cif_obj.write_file(
            str(filepath)
            )

    def block_index(self, block):

        blocks = [b for b in self.cif_obj]

        i = blocks.index(block)

        return i

    def iter_columns_from_suffix(self,
        column_suffix,
        block_name = None,
        ):

        for i_block, block in enumerate(self.cif_obj):

            if (block_name is not None) and (block.name != block_name):
                continue

            for i_item, item in enumerate(block):

                if item.loop is not None:

                    for i_col, col_name in enumerate(item.loop.tags):

                        if col_name.endswith(column_suffix):

                            col = block.find_loop(col_name)

                            yield (
                                block,
                                item,
                                col,
                                )

    def get_unique_column_values(self,
        block_name,
        column_name,
        ):

        block = self.cif_obj.find_block(
            block_name
            )

        assert block is not None

        values = []

        for item in block:

            if (item.loop is not None) and (column_name in item.loop.tags):

                table = block.item_as_table(item)

                col = table.find_column(column_name)

                values.extend(
                    list(col)
                    )

        return sorted(set(values))

    def prune_empty(self):

        # Prune empty blocks and items

        pass

    def replace_matching_values_by_column_name(self,
        column_name,
        old_value,
        new_value,
        block_name = None,
        ):

        for i_block, block in enumerate(self.cif_obj):

            if (block_name is not None) and (block.name != block_name):
                continue

            col = block.find_loop(column_name)

            if col is not None:

                for i, v in enumerate(col):

                    if v == old_value:

                        col[i] = new_value

    def replace_matching_values_by_column_suffix(self,
        column_suffix,
        old_value,
        new_value,
        block_name = None,
        ):

        for (block, item, col) in self.iter_columns_from_suffix(
            column_suffix = column_suffix,
            block_name = block_name,
            ):

            for i, v in enumerate(col):

                if v == old_value:

                    col[i] = new_value

    def remove_block(self,
        block_name,
        ):

        block = self.cif_obj.find_block(
            block_name,
            )

        i_block = self.block_index(block)

        del self.cif_obj[i_block]

    def remove_loop_row_by_value(self,
        block_name,
        column_name,
        row_value,
        ):
        # make iterative to remove multiple rows?

        block = self.cif_obj.find_block(
            block_name,
            )

        item = block.find_loop_item(
            column_name,
            )

        table = block.item_as_table(
            item,
            )

        col = table.find_column(
            column_name,
            )

        i_row = list(col).index(row_value)

        table.remove_row(i_row)

    def print_references_to(self, column=None, value=None):
        """
        Convenience function to find where in a document a column/value is
        """

        for block in self.cif_obj:

            for item in block:

                if (
                    column is not None
                    ) and (
                    column in item.loop.tags
                    ):

                    print 'column {} is in block {}'.format(
                        column,
                        block.name,
                        )

                if (
                    value is not None
                    ) and (
                    value in item.loop.values
                    ):

                    t = block.item_as_table(item)

                    for r in t:

                        vals = list(r)

                        if value in vals:

                            i = vals.index(value)

                            print 'value {} is in column {} in block {}'.format(
                                value,
                                t.column(i).tag,
                                block.name,
                                )

    # mods

    def all_mod_ids(self):

        return self.get_unique_column_values(
            block_name = self.MOD_LIST_NAME,
            column_name = (
                self.CHEM_MOD_ROOT + '.id'
                ),
            )

    # links

    def all_link_ids(self):

        return self.get_unique_column_values(
            block_name = self.LINK_LIST_NAME,
            column_name = (
                self.CHEM_LINK_ROOT + '.id'
                ),
            )

    def get_link_id_for_linking_atoms(self,
        atom_1_resname,
        atom_1_name,
        atom_2_resname,
        atom_2_name,
        ):

        ref_tuple = sorted([
            (atom_1_resname.strip(), atom_1_name.strip()),
            (atom_2_resname.strip(), atom_2_name.strip()),
            ])

        for (link_id, atom_1_tuple, atom_2_tuple) in self.iter_links_and_linking_atoms():

            if sorted([atom_1_tuple,atom_2_tuple]) == ref_tuple:

                return link_id

        return None

    def iter_links_and_linking_atoms(self):

        for link_id in self.all_link_ids():

            atom_1_tuple, atom_2_tuple = self.get_linking_atoms(link_id)

            yield (link_id, atom_1_tuple, atom_2_tuple)

    def get_linking_atoms(self, link_id):

        # get compound id refs as tuple of (res1, res2)

        resnames_tuple = None

        list_block = self.cif_obj.find_block(
            self.LINK_LIST_NAME
            )

        assert list_block is not None, 'Link block not found'

        list_item = list_block.find_loop_item('_chem_link.id')

        assert list_item is not None, 'Link list not found'

        i0 = list_item.loop.tags.index('_chem_link.id')
        i1 = list_item.loop.tags.index('_chem_link.comp_id_1')
        i2 = list_item.loop.tags.index('_chem_link.comp_id_2')

        list_table = list_block.item_as_table(list_item)

        for r in list_table:

            if r[i0] == link_id:

                resname_1 = r[i1]
                resname_2 = r[i2]

                resnames_tuple = (resname_1, resname_2)

                break

        if resnames_tuple is None:
            raise Exception('LinkID not found in link list')

        # now find the linking atoms for the link

        link_block = self.cif_obj.find_block(
            self.LINK_BLOCK_PRE + link_id
            )

        assert link_block is not None, 'Link description not found'

        link_item = link_block.find_loop_item('_chem_link_bond.link_id')

        assert link_item is not None, 'Link bonds not found'

        i0 = link_item.loop.tags.index('_chem_link_bond.link_id')
        i1 = link_item.loop.tags.index('_chem_link_bond.atom_1_comp_id')
        i2 = link_item.loop.tags.index('_chem_link_bond.atom_id_1')
        i3 = link_item.loop.tags.index('_chem_link_bond.atom_2_comp_id')
        i4 = link_item.loop.tags.index('_chem_link_bond.atom_id_2')

        link_table = link_block.item_as_table(link_item)

        for r in link_table:

            if r[i0] == link_id:

                i_resname_1 = int(r[i1]) - 1
                atom_name_1 = r[i2]
                i_resname_2 = int(r[i3]) - 1
                atom_name_2 = r[i4]

                resname_1 = resnames_tuple[i_resname_1]
                resname_2 = resnames_tuple[i_resname_2]

                return (
                    (resname_1.strip(), atom_name_1.strip()),
                    (resname_2.strip(), atom_name_2.strip()),
                    )

    # renaming

    def rename_block(self, old_name, new_name):

        block = self.cif_obj.find_block(old_name)

        assert block is not None, 'block not found: {}'.format(old_name)

        block.name = new_name

    def rename_mod(self, old_name, new_name):

        self.replace_matching_values_by_column_name(
            block_name = self.MOD_LIST_NAME,
            column_name = self.CHEM_MOD_ROOT + '.id',
            old_value = old_name,
            new_value = new_name,
            )

        for col_suffix in ['.mod_id_1', '.mod_id_2']:

            self.replace_matching_values_by_column_suffix(
                column_suffix = col_suffix,
                old_value = old_name,
                new_value = new_name,
                block_name = self.LINK_LIST_NAME, # Links -> mods
                )

        # now change the mod block itself

        old_block_name = self.MOD_BLOCK_PRE + old_name
        new_block_name = self.MOD_BLOCK_PRE + new_name

        self.rename_block(
            old_name = old_block_name,
            new_name = new_block_name,
            )

        for col_suffix in ['.mod_id']:

            self.replace_matching_values_by_column_suffix(
                column_suffix = col_suffix,
                old_value = old_name,
                new_value = new_name,
                block_name = new_block_name,
                )

    def rename_link(self, old_name, new_name):

        self.replace_matching_values_by_column_name(
            block_name = self.LINK_LIST_NAME,
            column_name = self.CHEM_LINK_ROOT + '.id',
            old_value = old_name,
            new_value = new_name,
            )

        # now change the link block itself

        old_block_name = self.LINK_BLOCK_PRE + old_name
        new_block_name = self.LINK_BLOCK_PRE + new_name

        self.rename_block(
            old_name = old_block_name,
            new_name = new_block_name,
            )

        for col_suffix in ['.link_id']:

            self.replace_matching_values_by_column_suffix(
                column_suffix = col_suffix,
                old_value = old_name,
                new_value = new_name,
                block_name = new_block_name,
                )


class CifIDGenerator:

    def __init__(self, seed=124, max_label_length=8):

        self.set_seed(seed)

        self.mod_regex = re.compile('[A-Z0-9]{3}m[0-9]+')

        self.link_regex = re.compile('[A-Z0-9]{3}-[A-Z0-9]{3}')

        self.link_regex_2 = re.compile('[A-Z0-9]{2}-[A-Z0-9]{2}-[0-9]{2}')

    def set_seed(self, seed):

        np.random.seed(seed)

    def __call__(self, name):

        ###

        if len(name) > 8:

            return self.random_label()

        ###

        match = self.mod_regex.match(name)

        if (match is not None) and (match.group() == name):

            return self.increment_mod(name)

        ###

        match = self.link_regex.match(name)

        if (match is not None) and (match.group() == name):

            return self.increment_link(name)

        ###

        match = self.link_regex_2.match(name)

        if (match is not None) and (match.group() == name):

            return self.increment_link_2(name)

        return self.random_label()

    def random_label(self):

        return ''.join(
            np.random.choice('ABCDEFGHIJKLMNOPQRSTUVWXYZ', size=8, replace=True)
            )

    def increment_mod(self, name):

        prefix = name[0:4]

        assert prefix[3] == 'm'

        i = int(name[4:])

        new_lab = prefix + str(i+1)

        return new_lab

    def increment_link(self, name):

        first, second = name.split('-')

        assert len(first) == 3
        assert len(second) == 3

        new_lab = (
            first[0:2] +
            '-' +
            second[0:2] +
            '-02' # because the original is the first
            )

        return new_lab

    def increment_link_2(self, name):

        first, second, i = name.split('-')

        assert len(first) == 2
        assert len(second) == 2

        i = int(i)

        new_lab = (
            first +
            '-' +
            second +
            '-{:02d}'.format(i+1)
            )

        return new_lab


class CifMerger:


    def __init__(self):

        self.get_new_id = CifIDGenerator()

    def __call__(self, cif_managers):

        logger('Removing duplicate compounds entries')
        self.remove_duplicate_compounds(cif_managers)

        logger('Renaming duplicated link names')
        self.resolve_link_clashes(cif_managers)

        logger('Renaming duplicated modification names')
        self.resolve_mod_clashes(cif_managers)

        #self.check_for_conflicts(cif_managers)

        logger('Merging cif objects')
        merged_cif = self.merge(cif_managers)

        return merged_cif

    def remove_duplicate_compounds(self, cif_managers):

        reference_cif = cif_managers[0]

        all_comp_ids = reference_cif.get_unique_column_values(
            block_name = 'comp_list',
            column_name = '_chem_comp.id',
            )

        for other_cif in cif_managers[1:]:

            for other_comp_id in other_cif.get_unique_column_values(
                block_name = 'comp_list',
                column_name = '_chem_comp.id',
                ):

                if other_comp_id in all_comp_ids:

                    logger(
                        'Removing repeated compound {comp_id} (from source: {src})'.format(
                            comp_id = other_comp_id,
                            src = str(other_cif.cif_obj.source),
                            )
                        )

                    other_cif.remove_loop_row_by_value(
                        block_name = 'comp_list',
                        column_name = '_chem_comp.id',
                        row_value = other_comp_id,
                        )

                    other_cif.remove_block(
                        block_name = 'comp_' + str(other_comp_id),
                        )

    def resolve_mod_clashes(self, cif_managers):

        reference_cif = cif_managers[0]

        all_ids = reference_cif.all_mod_ids()

        for other_cif in cif_managers[1:]:

            for cur_id in other_cif.all_mod_ids():

                if (cur_id in all_ids):

                    new_id = self.get_new_id(cur_id)
                    while new_id in all_ids:
                        new_id = self.get_new_id(new_id)

                    logger(
                        'Renaming repeated modification {cur_id} to {new_id} (from source: {src})'.format(
                            cur_id = cur_id,
                            new_id = new_id,
                            src = str(other_cif.cif_obj.source),
                            )
                        )

                    other_cif.rename_mod(
                        old_name = cur_id,
                        new_name = new_id,
                        )

                    all_ids.append(new_id)

                else:

                    all_ids.append(cur_id)

        return

    def resolve_link_clashes(self, cif_managers):

        reference_cif = cif_managers[0]

        all_ids = reference_cif.all_link_ids()

        for other_cif in cif_managers[1:]:

            for cur_id in other_cif.all_link_ids():

                if (cur_id in all_ids):

                    new_id = self.get_new_id(cur_id)
                    while new_id in all_ids:
                        new_id = self.get_new_id(new_id)

                    logger(
                        'Renaming repeated link {cur_id} to {new_id} (from source: {src})'.format(
                            cur_id = cur_id,
                            new_id = new_id,
                            src = str(other_cif.cif_obj.source),
                            )
                        )

                    other_cif.rename_link(
                        old_name = cur_id,
                        new_name = new_id,
                        )

                    all_ids.append(new_id)

                else:

                    all_ids.append(cur_id)

        return

    def merge(self, cif_managers):

        main_cif = cif_managers[0]

        for i_cif, other_cif in enumerate(cif_managers[1:]):

            logger(
                'Merging cif {n} (from source: {src}) into merged cif object'.format(
                    n = i_cif + 2, # (+1 for index to num, +1 for skipped first item)
                    src = str(other_cif.cif_obj.source),
                    )
                )

            for other_block in other_cif.cif_obj:

                main_block = main_cif.cif_obj.find_block(
                    other_block.name,
                    )

                if main_block is None:

                    logger(
                        'Copying block {block_id} to merged cif object'.format(
                            block_id = other_block.name,
                            )
                        )

                    main_cif.cif_obj.add_copied_block(
                        other_block
                        )

                else:

                    logger(
                        'Merging block {block_id} into existing block in merged cif object'.format(
                            block_id = other_block.name,
                            )
                        )

                    self.merge_blocks(
                        main_block = main_block,
                        other_block = other_block,
                        )

        return main_cif

    @staticmethod
    def merge_blocks(main_block, other_block):

        for other_item in other_block:

            other_loop = other_item.loop
            other_pair = other_item.pair

            if other_loop is not None:

                # Get matching column
                main_loop_col = main_block.find_loop(
                    other_loop.tags[0]
                    )

                if main_loop_col is None:

                    # new loop
                    main_block.add_item(
                        other_item
                        )
                else:

                    # existing loop
                    main_loop = main_loop_col.get_loop()

                    assert main_loop.tags == other_loop.tags, 'incompatible loops'

                    for other_row in other_block.item_as_table(other_item):

                        main_loop.add_row(
                            list(other_row)
                            )

            if other_pair is not None:

                # TODO
                raise Exception('not implemented')


#class LinkRecord:
#
#    def __init__(self, link_record):
#
#        self._parse(link_record)
#
#    def __str__(self):
#
#        return self.format_pdb()
#
#    def _parse(self, link_record):
#
#        if len(link_record) < 80:
#            link_record = (
#                link_record + ' '*(80-len(link_record))
#                )
#
#        self.link_type = (
#            link_record[0:6].strip(' ')
#            )
#
#        self.atom_1_name = (
#            link_record[12:16].strip(' ')
#            )
#        self.atom_1_altloc = (
#            link_record[16].strip(' ')
#            )
#        self.atom_1_resname = (
#            link_record[17:20].strip(' ')
#            )
#        self.atom_1_chain = (
#            link_record[21].strip(' ')
#            )
#        self.atom_1_resseq = (
#            link_record[22:26].strip(' ')
#            )
#        self.atom_1_inscode = (
#            link_record[26].strip(' ')
#            )
#
#        self.atom_2_name = (
#            link_record[42:46].strip(' ')
#            )
#        self.atom_2_altloc = (
#            link_record[46].strip(' ')
#            )
#        self.atom_2_resname = (
#            link_record[47:50].strip(' ')
#            )
#        self.atom_2_chain = (
#            link_record[51].strip(' ')
#            )
#        self.atom_2_resseq = (
#            link_record[52:56].strip(' ')
#            )
#        self.atom_2_inscode = (
#            link_record[56].strip(' ')
#            )
#
#        self.sym_op_1 = (
#            link_record[59:65].strip(' ')
#            )
#        self.sym_op_2 = (
#            link_record[66:72].strip(' ')
#            )
#
#        self.set_link_id(
#            link_record[72:80].strip(' ')
#            )
#
#    def format_pdb(self):
#
#        s = ''.join([
#            '{:<6}'.format(self.link_type),
#            ' '*6,
#            '{:^4}'.format(self.atom_1_name),
#            '{:>1}'.format(self.atom_1_altloc),
#            '{:>3}'.format(self.atom_1_resname),
#            '{:>2}'.format(self.atom_1_chain),
#            '{:>4}'.format(self.atom_1_resseq),
#            '{:>1}'.format(self.atom_1_inscode),
#            ' '*15,
#            '{:^4}'.format(self.atom_2_name),
#            '{:>1}'.format(self.atom_2_altloc),
#            '{:>3}'.format(self.atom_2_resname),
#            '{:>2}'.format(self.atom_2_chain),
#            '{:>4}'.format(self.atom_2_resseq),
#            '{:>1}'.format(self.atom_2_inscode),
#            ' '*2,
#            '{:<6}'.format(self.sym_op_1),
#            ' '*1,
#            '{:<6}'.format(self.sym_op_2),
#            '{:<8}'.format(self.link_id),
#            ])
#
#        assert len(s) == 80
#
#        return s
#
#    def set_link_id(self, link_id):
#
#        self.link_id = link_id.strip(' ')
#
#        if (self.link_id is not None) and self.link_id.strip(' '):
#            self.link_type = 'LINKR'
#        else:
#            self.link_type = 'LINK'
#
#    def set_link_id_from_cif(self, cif_manager):
#
#        link_id = cif_manager.get_link_id_for_linking_atoms(
#            self.atom_1_resname,
#            self.atom_1_name,
#            self.atom_2_resname,
#            self.atom_2_name,
#            )
#
#        if link_id is None:
#            raise Exception('No matching link found in cif')
#
#        self.set_link_id(link_id)
#
#
#class UpdateLinksFromCif:
#
#    LinkClass = LinkRecord
#
#    def __init__(self):
#
#        pass
#
#    def __call__(self,
#        cif_manager,
#        model,
#        ):
#
#        link_records = [
#            self.LinkClass(l)
#            for l in (
#                list(model.input.connectivity_annotation_section()) +
#                list(model.input.unknown_section())
#                )
#            if l.startswith('LINKR ') or l.startswith('LINK  ')
#            ]
#
#        for l in link_records:
#
#            l.set_link_id_from_cif(cif_manager)
#
#        return link_records
#
#
