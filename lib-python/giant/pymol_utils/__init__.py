import os, copy
import numpy

from itertools import cycle

from giant.pymol_utils.shapes import *


class PymolColourPalettes(object):


    set1 = ['carbon','cyan','lightmagenta','yellow','salmon','hydrogen','slate','orange']
    set2 = ['lime','deepteal','hotpink','yelloworange','violetpurple','grey70','marine','olive']
    set3 = ['smudge','teal','dirtyviolet','wheat','deepsalmon','lightpink','aquamarine','paleyellow']
    set4 = ['limegreen','skyblue','warmpink','limon','violet','bluewhite','greencyan','sand']
    set5 = ['forest','lightteal','darksalmon','splitpea','raspberry','grey50','deepblue','brown']


class _PymolScript(object):


    _colours = PymolColourPalettes.set1

    _styles = ['lines', 'cartoon', 'sticks', 'ribbon', 'dots', 'spheres', 'surface', 'mesh', 'ellipsoids']

    defaults = {'ray':     {'width':1024, 'height':786},
                'zoom':    {'buffer':1.0, 'complete':0},
                'rebuild': {'selection':'all','representation':'everything'}
               }

    def __init__(self, pretty_models=True, pretty_maps=True):

        # Previous object manipulated
        self.prev = None

        self.lines = []
        self._append(self._import_pml)
        self._append(self._import_cgo)

        if pretty_models:
            self.pretty_models()
        if pretty_maps:
            self.pretty_maps()

        self.colour_palette(self._colours)

    @classmethod
    def _build(cls):
        """Dynamically build the functions for this class"""

        for prefix in ['_cmd_']:
            for attribute_name in dir(cls):
                if not attribute_name.startswith(prefix):
                    continue
                func_name = attribute_name[len(prefix):]
                if hasattr(cls, func_name):
                    continue
                func = cls._make_method(func_name=func_name, template_name=attribute_name)
                setattr(cls, func_name, func)

    @classmethod
    def _make_method(cls, func_name, template_name):
        """Build custom functions for template population"""

        def _func(self, *args, **kw_args):
            defaults = copy.copy(self.defaults.get(func_name, {}))
            defaults.update(kw_args)
            template = getattr(self, template_name)
            cmd = self._format(template=template, *args, **defaults)
            self._append(cmd)
            return cmd

        return _func

    # ----------------------------------------------------------------------- #
    # Access/Modify Object Properties
    # ----------------------------------------------------------------------- #

    def colour_palette(self, colours):
        self._colours = colours
        self._colour_cycle = cycle(self._colours)

    def get_colour_palette(self):
        return self._colours

    def get_colour(self):
        return next(self._colour_cycle)

    # ----------------------------------------------------------------------- #
    # Append to/write script
    # ----------------------------------------------------------------------- #

    def _format(self, template, *args, **kw_args):
        assert kw_args.get('style') in [None,'everything']+self._styles
        return template.format(*args, **kw_args)

    def _append(self, thing):
        if isinstance(thing, list):     [self._append(t) for t in thing]
        elif isinstance(thing, str):    self.lines.append(thing)
        else: raise Exception('Invalid type: {}'.format(type(thing)))

    def append(self, thing):
        self._append(thing)

    def write_script(self, f_name, overwrite=False):
        if (overwrite is True) and os.path.exists(f_name):
            os.remove(f_name)
        assert not os.path.exists(f_name)
        assert f_name.endswith(self._file_type)
        with open(f_name, 'w') as fh:
            fh.write('\n'.join(self.lines)+'\n')
        return f_name

    # ----------------------------------------------------------------------- #
    # Modify the settings in PyMOL session
    # ----------------------------------------------------------------------- #

    def set(self, setting, *args, **kwargs):
        self._append(self._man_set \
                     .format(setting = repr(setting),
                             args    = ','.join([repr(a) for a in args]),
                             kwargs  = ','.join(['{}={}'.format(k, repr(kwargs[k])) for k in kwargs.keys()])) \
                     .replace(',,',',') \
                     .replace('(,','(') \
                     .replace(',)',')') \
                    )

    def set_normalise_maps(self, value=True):
        self.set("normalize_ccp4_maps", int(value))

    def custom(self, function, *args, **kwargs):
        self._append(self._man_custom \
                     .format(function = function,
                             args    = ','.join([repr(a) for a in args]),
                             kwargs  = ','.join(['{}={}'.format(k, repr(kwargs[k])) for k in kwargs.keys()])) \
                     .replace(',,',',') \
                     .replace('(,','(') \
                     .replace(',)',')') \
                    )

    # ----------------------------------------------------------------------- #
    # Load and fetch structures and maps
    # ----------------------------------------------------------------------- #

    def fetch_pdb(self, pdb_code):
        pass

    def fetch_map(self, pdb_code, map_type='2fofc'):
        pass

    def load_pdb(self, f_name, obj=None, state=0):
        """Load pdb file f_name"""
        if obj is None: obj = os.path.basename(os.path.splitext(f_name)[0])
        self._append(self._load_basic.format(f_name=f_name, obj=obj, state=state))
        return obj

    def load_map(self, f_name, obj=None, state=0):
        if obj is None: obj = os.path.basename(os.path.splitext(f_name)[0])
        self._append(self._load_basic.format(f_name=f_name, obj=obj, state=state))
        return obj

    def make_mesh(self, obj, mesh_suffix='.mesh', contour_level=1, colour=None):
        mesh_name = obj+mesh_suffix
        self._append(self._isomesh.format(obj=mesh_name, map=obj, contour_level=contour_level))
        if colour is not None: self.colour(obj=mesh_name, colour=colour)
        return mesh_name

    # ----------------------------------------------------------------------- #
    # Quick addition of Shape objects
    # ----------------------------------------------------------------------- #

    def add_generic_object(self, item, obj, state=0):
        self._append(item.as_cmd(obj=obj, state=state))

    def add_shape(self, shape, obj, state=0):
        self.add_generic_object(item=shape, obj=obj, state=state)

    def add_shapes(self, cgo_shapes, obj, state=0):
        cgo_str = ''
        for c in cgo_shapes:
            c_str = c.as_string().strip('[]')
            cgo_str += c_str + ', '
        cgo_str = '[' + cgo_str + ']'
        cmd = 'cmd.load_cgo(%s, "%s", state=%s)' % (cgo_str, obj, state)
        self._append(cmd)

    # ----------------------------------------------------------------------- #
    # Multi-command defaults
    # ----------------------------------------------------------------------- #

    def pretty_models(self):
        self.set("cartoon_side_chain_helper", 1)

    def pretty_maps(self):
        pass

    # ----------------------------------------------------------------------- #
    # Visualisation
    # ----------------------------------------------------------------------- #

    def repr_as(self, *args, **kw_args):
        return self.show_as(*args, **kw_args)

    def repr_show(self,*args, **kw_args):
        return self.show(*args, **kw_args)

    def repr_hide(self, style='everything', *args, **kw_args):
        kw_args['style'] = style
        return self.hide(*args, **kw_args)

    def colour(self, obj, colour=None):
        if colour is None: colour = self.get_colour()
        cmd = self._format(template=self._man_colour, colour=repr(colour), obj=repr(obj))
        self._append(cmd)
        return colour

    def colour_by_element(self, obj, carbon_colour='green'):
        cmd = self._format(template=self._util_colour_by[carbon_colour], obj=obj)
        self._append(cmd)
        return obj

    def chainbow(self, obj):
        cmd = self._format(template=self._util_chainbow, obj=obj)
        self._append(cmd)
        return obj

    # ----------------------------------------------------------------------- #

class PythonScript(_PymolScript):


    _file_type = '.py'

    _import_pml = 'from pymol import cmd, util'
    _import_cgo = 'from pymol.cgo import *'

    _load_basic = 'cmd.load("{f_name}","{obj}",state={state})'
    _isomesh    = 'cmd.isomesh("{obj}", "{map}", {contour_level})'

    _man_set        = 'cmd.set({setting},{args},{kwargs})'
    _man_colour     = 'cmd.color({colour}, selection={obj})'
    _man_custom     = 'cmd.{function}({args},{kwargs})'

    _cmd_quit       = 'cmd.quit()'
    _cmd_ray        = 'cmd.ray({width},{height})'
    _cmd_select     = 'cmd.select("{obj}", "{selection}")'
    _cmd_show_as    = 'cmd.show_as("{style}", "{obj}")'
    _cmd_show       = 'cmd.show("{style}", "{obj}")'
    _cmd_hide       = 'cmd.hide("{style}", "{obj}")'
    _cmd_orient     = 'cmd.orient("{obj}")'
    _cmd_zoom       = 'cmd.zoom("{obj}", buffer={buffer}, complete={complete})'
    _cmd_png        = 'cmd.png("{f_name}")'
    _cmd_disable    = 'cmd.disable("{obj}")'
    _cmd_enable     = 'cmd.enable("{obj}")'

    _cmd_label      = 'cmd.label("{selection}", "{expression}")'
    _cmd_rebuild    = 'cmd.rebuild(selection="{selection}", representation="{representation}")'

    _cmd_ramp_new   = 'cmd.ramp_new("{obj}", map_name="{mol_or_map_obj}", range={range}, color={color})'

    _util_colour_by = {'black'    : 'util.cbab("{obj}")',
                       'green'    : 'util.cbag("{obj}")',
                       'magenta'  : 'util.cbam("{obj}")',
                       'purple'   : 'util.cbap("{obj}")',
                       'white'    : 'util.cbaw("{obj}")',
                       'cyan'     : 'util.cbac("{obj}")',
                       'pink'     : 'util.cbak("{obj}")',
                       'orange'   : 'util.cbao("{obj}")',
                       'salmon'   : 'util.cbas("{obj}")',
                       'yellow'   : 'util.cbay("{obj}")'}

    _util_chainbow = 'util.chainbow("{obj}")'

    def change_into_directory(self, path):
        self._append('import os')
        self._append('os.chdir({})'.format(repr(path)))

    def change_into_directory_maybe(self, path):
        self._append('import os')
        self._append('if os.path.exists({d}):\n    os.chdir({d})'.format(d=repr(path)))

    @classmethod
    def run(cls, script):
        assert script.endswith(cls._file_type)
        from giant.dispatcher import Dispatcher
        prog = Dispatcher('pymol')
        prog.extend_args(['-k', '-q', '-c', '-Q', '-s', script[:-len(cls._file_type)]+'.log', '-r', script])
        prog.run()
        return prog

# Make the custom functions for the class at runtime
PythonScript._build()
