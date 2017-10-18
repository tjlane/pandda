
from iotbx.pdb.hierarchy import input as ih

from bamboo.pymol_utils import PythonScript
from giant.structure.select import get_select_function
from giant.structure.formatting import PymolSelection, ShortLabeller

def auto_chain_images(structure_filename, output_prefix, selection='protein', style='cartoon', **kw_args):
    filter_func = get_select_function(selection)
    return selection_images(structure_filename  = structure_filename,
                            output_prefix       = output_prefix+'-chain_',
                            selections  = [PymolSelection.format(c) for c in filter_func(ih(structure_filename).hierarchy).chains()],
                            labels      = [ShortLabeller.format(c)  for c in filter_func(ih(structure_filename).hierarchy).chains()],
                            style = style, **kw_args)

def auto_residue_images(structure_filename, output_prefix, selection='protein', style='sticks', **kw_args):
    filter_func = get_select_function(selection)
    return selection_images(structure_filename  = structure_filename,
                            output_prefix       = output_prefix+'-residue_',
                            selections  = [PymolSelection.format(r) for r in filter_func(ih(structure_filename).hierarchy).residue_groups()],
                            labels      = [ShortLabeller.format(r)  for r in filter_func(ih(structure_filename).hierarchy).residue_groups()],
                            style = style, **kw_args)

def selection_images(structure_filename,
                     output_prefix,
                     selections,
                     labels = None,
                     style  = 'sticks',
                     ray_trace      = True,
                     run_script     = True,
                     delete_script  = True,
                    ):

    if labels is None:
        labels = map(str,range(1, len(selections)+1))

    # Create script object
    s = PythonScript(pretty_models=False, pretty_maps=False)
    s.set('opaque_background', 'off')
    s.set('orthoscopic', 1)
    # Read in structure and hide all atoms
    s.load_pdb(f_name=structure_filename, obj='input_structure')

    # Set the styles
    style = style.split('+')
    # Set colors
    s.colour(obj='all', colour="grey50")
    s.show_as(obj='all', style='lines')

    png_filenames = []

    for i, selection in enumerate(selections):
        s.select(obj='sele', selection=selection)
        s.orient(obj='sele')
        s.zoom(obj='sele', buffer=1, complete=0)
        s.colour_by_element(obj='sele', carbon_colour='green')
        for sty in style:
            s.show(obj='sele', style=sty)
        if ray_trace:
            s.ray()
        png_name = s.png(f_name=output_prefix+labels[i]+'.png')
        png_filenames.append(png_name)
        s.colour(obj='sele', colour="grey50")
        s.show_as(obj='sele', style='lines')

    f_name = s.write_script(output_prefix+'.py')

    if run_script is True:
        s.run(f_name)

    return png_filenames
