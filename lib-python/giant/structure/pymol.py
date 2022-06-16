import giant.logs as lg
logger = lg.getLogger(__name__)

import os, itertools, collections

from iotbx.pdb.hierarchy import input as ih

from giant.pymol_utils import PythonScript
from giant.structure.select import get_select_function
from giant.structure.formatting import PymolSelection, ShortLabeller

def _load_and_filter(structure_filename, selection):
    s_h = ih(structure_filename).hierarchy
    c = s_h.atom_selection_cache()
    s = c.selection(selection)
    return s_h.select(s)

def auto_chain_images(
        structure_filename,
        output_prefix,
        selection = 'not water',
        style = 'cartoon',
        het_style = None,
        **kw_args
        ):
    h = _load_and_filter(structure_filename, selection)
    output_files = selection_images(
        structure_filename = structure_filename,
        output_prefix = output_prefix+'-chain_',
        selections = [PymolSelection.format(c) for c in h.chains()],
        labels     = [ShortLabeller.format(c)  for c in h.chains()],
        style = style, het_style=het_style, **kw_args)
    return collections.OrderedDict(zip([c.id for c in h.chains()], output_files))

def auto_residue_images(
        structure_filename,
        output_prefix,
        selection = 'not water',
        style = 'sticks',
        het_style = None,
        **kw_args
        ):
    h = _load_and_filter(structure_filename, selection)
    output_files = selection_images(
        structure_filename = structure_filename,
        output_prefix = output_prefix+'-residue_',
        selections = [PymolSelection.format(r) for r in h.residue_groups()],
        labels     = [ShortLabeller.format(r)  for r in h.residue_groups()],
        style = style, het_style=het_style, **kw_args)
    return collections.OrderedDict(zip([r.id_str() for r in h.residue_groups()], output_files))

def selection_images(
        structure_filename,
        output_prefix,
        selections,
        labels = None,
        style = 'sticks',
        colours = None,
        hide_rest = True,
        ray_trace = True,
        settings = [],
        run_script = True,
        delete_script = False,
        width = 800,
        height = 600,
        colour_selections = None,
        het_style = None,
        ):

    if labels is None:
        labels = list(map(str,range(1, len(selections)+1)))

    assert (colours is None) or (colours in ['bfactor','chainbow']) or isinstance(colours, list)
    if (colours is None):
        colours = ['green']
    if isinstance(colours, list):
        colours = itertools.cycle(colours)

    if colour_selections is not None:
        assert colours != 'bfactor'
        assert colours != 'chainbow'

    # Create script object
    s = PythonScript(pretty_models=False, pretty_maps=False)
    s.custom('bg_color', 'white')
    s.set('opaque_background', 0)
    s.set('ray_opaque_background', 0)
    s.set('orthoscopic', 1)
    # Apply custom global commands
    for cmd in settings:
        s.set(*cmd)
    # Read in structure and hide all atoms
    s.load_pdb(f_name=structure_filename, obj='input_structure')
    # Set styles
    style = style.split('+')
    if het_style:
        het_style = het_style.split('+')
    else:
        het_style = []

    png_filenames = []

    for i, selection in enumerate(selections):
        # Reset the structure to all look the same
        s.show_as(obj='all', style=style[0])
        s.colour(obj='all', colour="grey90")
        for sty in het_style:
            s.show(obj='het', style=sty)
        # Apply custom views to the selection
        s.select(obj='sele', selection=selection)
        for sty in style:
            s.show(obj='sele', style=sty)
        s.orient(obj='sele')
        s.zoom(obj='sele', buffer=0.0, complete=1)
        # Custom colouring
        if colour_selections:
            for sel in colour_selections:
                s.colour(obj="sele and ({})".format(sel), colour=next(colours))
        elif colours == 'bfactor':
            s.custom('spectrum', expression='b', palette='blue_red', selection='sele and (b>-1)')
        elif colours == 'chainbow':
            s.chainbow("sele")
        else:
            s.colour_by_element(obj='sele', carbon_colour=next(colours))
        # Hide
        if hide_rest:
            s.hide(obj='not sele', style='everything')
        # Ray trace
        if ray_trace:
            s.ray(height=height, width=width)
        f_name = output_prefix+labels[i]+'.png'
        s.png(f_name=f_name)
        png_filenames.append(f_name)

    f_name = s.write_script(output_prefix+'.py')
    l_name = f_name.replace('.py', '.log')
    assert not os.path.exists(l_name)

    if (run_script is True):
        o = s.run(f_name)
        o.write_output(l_name)
        for f in png_filenames:
            if not os.path.exists(f):
                delete_script = False

    if (delete_script is True):
        if os.path.exists(f_name):
            os.remove(f_name)
        if os.path.exists(l_name):
            os.remove(l_name)

    return png_filenames
