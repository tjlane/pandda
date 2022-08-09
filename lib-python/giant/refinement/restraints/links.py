import giant.logs as lg
logger = lg.getLogger(__name__)


class LinkRecord:

    def __init__(self, link_record):

        self._parse(link_record)

    def __str__(self):

        return self.format_pdb()

    def _parse(self, link_record):

        if len(link_record) < 80:
            link_record = (
                link_record + ' '*(80-len(link_record))
                )

        self.link_type = (
            link_record[0:6].strip(' ')
            )

        self.atom_1_name = (
            link_record[12:16].strip(' ')
            )
        self.atom_1_altloc = (
            link_record[16].strip(' ')
            )
        self.atom_1_resname = (
            link_record[17:20].strip(' ')
            )
        self.atom_1_chain = (
            link_record[21].strip(' ')
            )
        self.atom_1_resseq = (
            link_record[22:26].strip(' ')
            )
        self.atom_1_inscode = (
            link_record[26].strip(' ')
            )

        self.atom_2_name = (
            link_record[42:46].strip(' ')
            )
        self.atom_2_altloc = (
            link_record[46].strip(' ')
            )
        self.atom_2_resname = (
            link_record[47:50].strip(' ')
            )
        self.atom_2_chain = (
            link_record[51].strip(' ')
            )
        self.atom_2_resseq = (
            link_record[52:56].strip(' ')
            )
        self.atom_2_inscode = (
            link_record[56].strip(' ')
            )

        self.sym_op_1 = (
            link_record[59:65].strip(' ')
            )
        self.sym_op_2 = (
            link_record[66:72].strip(' ')
            )

        self.set_link_id(
            link_record[72:80].strip(' ')
            )

    def format_pdb(self):

        s = ''.join([
            '{:<6}'.format(self.link_type),
            ' '*6,
            '{:^4}'.format(self.atom_1_name),
            '{:>1}'.format(self.atom_1_altloc),
            '{:>3}'.format(self.atom_1_resname),
            ' ',
            '{:>1}'.format(self.atom_1_chain),
            '{:>4}'.format(self.atom_1_resseq),
            '{:>1}'.format(self.atom_1_inscode),
            ' '*15,
            '{:^4}'.format(self.atom_2_name),
            '{:>1}'.format(self.atom_2_altloc),
            '{:>3}'.format(self.atom_2_resname),
            ' ',
            '{:>1}'.format(self.atom_2_chain),
            '{:>4}'.format(self.atom_2_resseq),
            '{:>1}'.format(self.atom_2_inscode),
            ' '*2,
            '{:>6}'.format(self.sym_op_1),
            ' '*1,
            '{:>6}'.format(self.sym_op_2),
            '{:>8}'.format(self.link_id),
            ])

        assert len(s) == 80

        return s

    def set_link_id(self, link_id):

        self.link_id = link_id.strip(' ')

        if (self.link_id is not None) and self.link_id.strip(' '):
            self.link_type = 'LINKR'
        else:
            self.link_type = 'LINK'

    def set_link_id_from_cif(self, cif_manager):

        link_id = cif_manager.get_link_id_for_linking_atoms(
            self.atom_1_resname,
            self.atom_1_name,
            self.atom_2_resname,
            self.atom_2_name,
            )

        if link_id is None:
            raise Exception('No matching link found in cif')

        self.set_link_id(link_id)


class UpdateLinksFromCif:

    LinkClass = LinkRecord

    def __init__(self):

        pass

    def __call__(self,
        cif_manager,
        model,
        ):

        link_records = [
            self.LinkClass(l)
            for l in (
                list(model.input.connectivity_annotation_section()) +
                list(model.input.unknown_section())
                )
            if l.startswith('LINKR ') or l.startswith('LINK  ')
            ]

        link_records = sorted(
            link_records,
            key = lambda l: (
                str(l.atom_1_chain),
                int(l.atom_1_resseq),
                str(l.atom_1_inscode),
                str(l.atom_1_altloc),
                str(l.atom_2_chain),
                int(l.atom_2_resseq),
                str(l.atom_2_inscode),
                str(l.atom_2_altloc),
                str(l.atom_1_name),
                str(l.atom_2_name),
                ),
            )

        for l in link_records:

            logger('Model Link:  {l}'.format(l=str(l)))

            l.set_link_id_from_cif(cif_manager)

            logger('  Updated -> {l}\n'.format(l=str(l)))

        return link_records


