import iotbx.pdb.hierarchy as iotbx_pdbh


class _Selection(object):


    _labels = ('model','chain','resseq','icode','resname','altloc','name')

    def __call__(self, obj):
        return self.format(obj)

    @classmethod
    def format(cls, obj):
        return cls.join_and([s for s in cls._format_any_to_list(obj) if s not in cls.remove])

    @classmethod
    def join_or(cls, strings, extra_join=''):
        return (cls._join_or+extra_join).join([cls._sep_or.format(s) for s in strings])
    @classmethod
    def join_and(cls, strings, extra_join=''):
        return (cls._join_and+extra_join).join([cls._sep_and.format(s) for s in strings])
    @classmethod
    def join_custom(cls, strings, join):
        if '\n' in join:
            return join.join([s.replace('\n', join) for s in strings])
        else:
            return join.join(strings)

    @classmethod
    def _format_any_to_list(cls, obj):
        """Return the relevant label for a supplied hierarchy/atom object"""
        if   isinstance(obj, iotbx_pdbh.model):
            s = cls._format_mo(obj)
        elif isinstance(obj, iotbx_pdbh.chain):
            s = cls._format_ch(obj)
        elif isinstance(obj, iotbx_pdbh.residue_group):
            s = cls._format_rg(obj)
        elif isinstance(obj, iotbx_pdbh.atom_group):
            s = cls._format_ag(obj)
        elif isinstance(obj, iotbx_pdbh.conformer):
            s = cls._format_co(obj)
        elif isinstance(obj, iotbx_pdbh.residue):
            s = cls._format_re(obj)
        elif isinstance(obj, iotbx_pdbh.atom):
            if hasattr(obj, 'chain_id'): s = cls._format_al(obj)
            else:                        s = cls._format_at(obj)
        elif isinstance(obj, iotbx_pdbh.atom_with_labels):
            s = cls._format_al(obj)
        elif isinstance(obj, dict):
            s = cls._format_dict(obj)
        else:
            raise Exception('Invalid object type provided: {}'.format(type(obj)))
        return s

    @classmethod
    def _format_mo(cls, obj):
        return [ cls.model.format(obj.id) ]
    @classmethod
    def _format_ch(cls, obj):
        return cls._format_any_to_list(obj.parent()) + \
                [ cls.chain.format(obj.id),
                ]
    @classmethod
    def _format_rg(cls, obj):
        return cls._format_any_to_list(obj.parent()) + \
                [ cls.resseq.format(obj.resseq),
                  cls.icode.format(obj.icode),
                ]
    @classmethod
    def _format_ag(cls, obj):
        return cls._format_any_to_list(obj.parent()) + \
                [ cls.resname.format(obj.resname),
                  cls.altloc.format(obj.altloc),
                ]
    @classmethod
    def _format_co(cls, obj):
        return cls._format_any_to_list(obj.parent()) + \
                [ cls.altloc.format(obj.altloc),
                ]
    @classmethod
    def _format_re(cls, obj):
        return cls._format_any_to_list(obj.parent()) + \
                [ cls.resname.format(obj.resname),
                  cls.resseq.format(obj.resseq),
                  cls.icode.format(obj.icode),
                ]
    @classmethod
    def _format_at(cls, obj):
        return cls._format_any_to_list(obj.parent()) + \
                [ cls.name.format(obj.name),
                ]
    @classmethod
    def _format_al(cls, obj):
        return [cls.model.format(obj.model_id),
                cls.chain.format(obj.chain_id),
                cls.resseq.format(obj.resseq),
                cls.icode.format(obj.icode),
                cls.resname.format(obj.resname),
                cls.altloc.format(obj.altloc),
                cls.name.format(obj.name),
               ]
    @classmethod
    def _format_dict(cls, obj_dict):
        return [ cls.__dict__.get(l).format(obj_dict.get(l)) for l in cls._labels if (l in obj_dict) ]


class Labeller(_Selection):


    _join_and = '-'
    _join_or  = None    # Refmac has no equivalent atom selection syntax
    _sep_and  = '{}'
    _sep_or   = None    # Refmac has no equivalent atom selection syntax

    remove = ['',' ']

    model       = 'mod({})'
    chain       = 'chn({})'
    resseq      = 'res({})'
    icode       = 'ins({})'
    resname     = '{}'
    altloc      = 'alt({})'
    name        = '[{}]'


class ShortLabeller(Labeller):


    model       = ''
    chain       = '{}'
    resseq      = '{}'
    icode       = '{}'
    resname     = '{}'
    altloc      = '({})'
    name        = '[{}]'


class GenericSelection(_Selection):


    _join_and = '-'
    _join_or  = None
    _sep_and  = '{}'
    _sep_or   = None

    remove = []

    model       = '{}'
    chain       = '{}'
    resseq      = '{}'
    icode       = '{}'
    resname     = '{}'
    altloc      = '{}'
    name        = '{}'

    @classmethod
    def to_str(cls, obj):
        obj_dict = cls.to_dict(obj)
        obj_list = ['{}({})'.format(l, obj_dict.get(l,'*')) for l in cls._labels]
        return ','.join(obj_list)

    @classmethod
    def to_dict(cls, obj):
        """Format any object to a dictionary (except atoms)"""

        if   isinstance(obj, iotbx_pdbh.model):
            labs = ('model')
            info = cls._format_mo(obj)
        elif isinstance(obj, iotbx_pdbh.chain):
            labs = ('model','chain')
            info = cls._format_ch(obj)
        elif isinstance(obj, iotbx_pdbh.residue_group):
            labs = ('model','chain','resseq','icode')
            info = cls._format_rg(obj)
        elif isinstance(obj, iotbx_pdbh.atom_group):
            labs = ('model','chain','resseq','icode','resname','altloc')
            info = cls._format_ag(obj)
        elif isinstance(obj, iotbx_pdbh.conformer):
            labs = ('model','chain','altloc')
            info = cls._format_co(obj)
        elif isinstance(obj, iotbx_pdbh.residue):
            labs = ('model','chain','altloc','resname','resseq','icode')
            info = cls._format_re(obj)
        elif isinstance(obj, iotbx_pdbh.atom):
            labs = ('model','chain','resseq','icode','resname','altloc','name')
            if hasattr(obj, 'chain_id'): info = cls._format_al(obj)
            else:                        info = cls._format_at(obj)
        elif isinstance(obj, iotbx_pdbh.atom_with_labels):
            labs = ('model','chain','resseq','icode','resname','altloc','name')
            info = cls._format_al(obj)
        else:
            raise Exception('Invalid object type provided: {}'.format(type(obj)))

        assert len(labs) == len(info)
        return dict(zip(labs, info))


class RefmacSelection(_Selection):


    _join_and = ' '
    _join_or  = None    # Refmac has no equivalent atom selection syntax
    _sep_and  = '{}'
    _sep_or   = None    # Refmac has no equivalent atom selection syntax

    remove = ['','model ','alte ', 'insc  ']

    model       = 'model {}'
    chain       = 'chain {}'
    resseq      = 'resi {}'
    icode       = 'insc {}'
    resname     = ''
    altloc      = 'alte {}'
    name        = 'atom {}'


class PhenixSelection(_Selection):


    _join_and = ' and '
    _join_or  = ' or '
    _sep_and  = '{}'
    _sep_or   = '({})'

    remove = ['','model ']

    model       = "model {}"
    chain       = "chain '{:1}'"
    resseq      = "resseq {}"
    icode       = "icode '{:1}'"
    resname     = "resname '{}'"
    altloc      = "altid '{:1}'"
    name        = "name {}"


class PymolSelection(_Selection):


    _join_and = ' and '
    _join_or  = ' or '
    _sep_and  = '{}'
    _sep_or   = '({})'

    remove = ['','model ']

    model       = "model {}"
    chain       = "chain {:1}"
    resseq      = "resi {}"
    icode       = ""
    resname     = "resn '{}'"
    altloc      = "alt '{}'"
    name        = "name {}"


labeller = Labeller()
short_labeller = ShortLabeller()
