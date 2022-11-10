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

        self.set_link_type(
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

        if self.link_type == "LINK": 
            self.set_link_dist(
                link_record[73:78].strip(' ')
                )    
        else: 
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
            (
                '{:>8}'.format(
                    self.link_id
                    )
                if 
                self.link_type == "LINKR"
                else 
                ' {:>5}  '.format(
                    self.link_dist if (self.link_dist is not None) else ''
                    )
                ),
            ])

        assert len(s) == 80

        return s

    def set_link_type(self, link_type):

        link_type = link_type.strip()

        assert link_type in ["LINK","LINKR"]

        self.link_type = link_type

    def set_link_dist(self, link_dist): 

        if (link_dist is None) or (link_dist == ''):
            self.link_dist = None
            return

        self.link_dist = float(link_dist)
        self.set_link_type("LINK")
        self.set_link_id(None)

    def set_link_id(self, link_id):

        if link_id is None: 
            self.link_id = None
            return

        self.link_id = link_id.strip(' ')
        self.set_link_type("LINKR")
        self.set_link_dist(None)

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

        for l in link_records:

            logger('\nModel Link:  {l}'.format(l=str(l)))

            try:
                l.set_link_id_from_cif(cif_manager)
                logger('  Updated -> {l}'.format(l=str(l)))
            except Exception as e: 
                w = "{}: {}".format(str(e),str(l))
                logger.warning(str(w))
                continue

        link_records = sorted(
            link_records,
            key = lambda l: (
                str(l.link_type),
                str(l.atom_1_chain),
                str(l.atom_2_chain),
                int(l.atom_1_resseq),
                str(l.atom_1_inscode),
                str(l.atom_1_altloc),
                int(l.atom_2_resseq),
                str(l.atom_2_inscode),
                str(l.atom_2_altloc),
                str(l.atom_1_name),
                str(l.atom_2_name),
                ),
            )

        return link_records


