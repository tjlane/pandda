import os, copy
import numpy

from itertools import cycle

from bamboo.common.command import CommandManager
from bamboo.pymol_utils.shapes import *


class PymolColourPalettes:


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
               }

    def __init__(self, pretty_models=True, pretty_maps=True):

        # Previous object manipulated
        self.prev = None

        self.lines = []
        self._append(self._import_pml)
        self._append(self._import_cgo)

        if pretty_models: self.pretty_models()
        if pretty_maps:   self.pretty_maps()

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
        return self._colour_cycle.next()

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
                     .format(setting = setting,
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

    def load_pdb(self, f_name, obj=None):
        """Load pdb file f_name"""
        if obj is None: obj = os.path.basename(os.path.splitext(f_name)[0])
        self._append(self._load_basic.format(f_name=f_name, obj=obj))
        return obj

    def load_map(self, f_name, obj=None):
        if obj is None: obj = os.path.basename(os.path.splitext(f_name)[0])
        self._append(self._load_basic.format(f_name=f_name, obj=obj))
        return obj

    def make_mesh(self, obj, mesh_suffix='.mesh', contour_level=1, colour=None):
        mesh_name = obj+mesh_suffix
        self._append(self._isomesh.format(obj=mesh_name, map=obj, contour_level=contour_level))
        if colour is not None: self.colour(obj=mesh_name, colour=colour)
        return mesh_name

    # ----------------------------------------------------------------------- #
    # Quick addition of Shape objects
    # ----------------------------------------------------------------------- #

    def add_shape(self, shape, obj):
        self._append(shape.as_cmd(name=obj))

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
        cmd = self._format(template=self._man_colour, colour=colour, obj=obj)
        self._append(cmd)
        return colour

    def colour_by_element(self, obj, carbon_colour='green'):
        cmd = self._format(template=self._util_colour_by[carbon_colour], obj=obj)
        self._append(cmd)
        return obj

    # ----------------------------------------------------------------------- #

class PythonScript(_PymolScript):


    _file_type = '.py'

    _import_pml = 'from pymol import cmd, util'
    _import_cgo = 'from pymol.cgo import *'

    _load_basic = 'cmd.load("{f_name}","{obj}")'
    _isomesh    = 'cmd.isomesh("{obj}", "{map}", {contour_level})'

    _man_set        = 'cmd.set("{setting}",{args},{kwargs})'
    _man_colour     = 'cmd.color("{colour}", "{obj}")'
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

    @classmethod
    def run(cls, script):
        assert script.endswith(cls._file_type)
        c = CommandManager('pymol')
        c.add_command_line_arguments(['-k', '-q', '-c', '-Q', '-s', script[:-len(cls._file_type)]+'.log', '-r', script])
        c.run()
        return c

# Make the custom functions for the class at runtime
PythonScript._build()
