import os, glob
import pandas

from scitbx.array_family import flex

from giant.structure.formatting import Labeller, ShortLabeller

from bamboo.html import png2base64str
from pandemic.html import PANDEMIC_HTML_ENV

def write_adp_summary(parameterisation, out_dir=None):
    """Write an overall, and a level-by-level summary of the parameterised ADPs"""
    if out_dir is None: out_dir = parameterisation.out_dir
    # Get the html template
    template = PANDEMIC_HTML_ENV.get_template('adp_summary.html')
    # Create output file path
    out_file = os.path.abspath(os.path.join(out_dir, 'results.html'))

    # ===========================================================>
    # Construct the data object to populate the template
    output_data = {}
    output_data['header'] = 'pandemic.adp summary'
    output_data['title'] = 'PanDEMIC ADP Parameterisation'
    output_data['introduction'] = 'ADP parameterisation summary from pandemic.adp'
    # ===========================================================>
    # Construct the tabs
    output_data['tabs'] = []

    # Extract the fitter
    fitter = parameterisation.fitter
    table = parameterisation.table

    ##########################################################################################
    # Create overview tab
    ##########################################################################################
    tab = {'id'         : 'leveloverview',
           'active'     : True,
           'type'       : 'overview',
           'short_name' : 'Overview',
           'long_name'  : 'Overview of Level-by-level TLS parameterisation',
           'description': 'Composed of {} levels.'.format(len(fitter.levels)),
           'rows'       : [],
          }
    output_data['tabs'].append(tab)
    # Format each TLS level
    for i_level, (level_num, level_lab, level) in enumerate(fitter):
        chain_image = os.path.join(parameterisation.out_dir, 'model/pymol/level_{}-all-modes-chain_A.png'.format(level_num))
        stack_image = os.path.join(parameterisation.out_dir, 'model/level_{}-all-modes-TLS-chain_A.png'.format(level_num))
        row = {'id'             : 'Level {} of {} ({})'.format(level_num, len(fitter.levels),level_lab),
               'text'           : '{} atoms.'.format('X'),
               'chain_image'    : 'data:image/png;base64,{}'.format(png2base64str(chain_image)) if os.path.exists(chain_image) else '',
               'stack_image'    : 'data:image/png;base64,{}'.format(png2base64str(stack_image)) if os.path.exists(stack_image) else '',
              }
        tab['rows'].append(row)
    # Format residual level
    chain_image = os.path.join(parameterisation.out_dir, 'model/pymol/residual-chain_A.png')
    stack_image = os.path.join(parameterisation.out_dir, 'model/all-residual-chain_A.png')
    row = {'id'             : 'Final Level (residual)',
           'text'           : '{} atoms.'.format('X'),
           'chain_image'    : 'data:image/png;base64,{}'.format(png2base64str(chain_image)) if os.path.exists(chain_image) else '',
           'stack_image'    : 'data:image/png;base64,{}'.format(png2base64str(stack_image)) if os.path.exists(stack_image) else '',
          }
    tab['rows'].append(row)

    ##########################################################################################
    # Create tab for each level
    ##########################################################################################
    for i_level, (level_num, level_lab, level) in enumerate(fitter):
        # Create dictionary for this tab and add to tab_list
        tab = {'id'         : 'lvl{}'.format(level_num),
               'type'       : 'level',
               'short_name' : 'Level {}'.format(level_num),
               'long_name'  : 'Level {} ({})'.format(level_num, level_lab),
               'description': 'Level {} of {}. '.format(level_num, len(fitter.levels))+\
                              'Composed of {} groups'.format(level.n_groups()),
               'rows'       : [],
              }
        output_data['tabs'].append(tab)
        # Read in the TLS models and amplitudes for this level
        tls_models     = pandas.read_csv(os.path.join(parameterisation.out_dir, 'model/csvs/tls_models_level_{:04d}.csv'.format(level_num))).set_index(['group','mode'])
        tls_amplitudes = pandas.read_csv(os.path.join(parameterisation.out_dir, 'model/csvs/tls_amplitudes_level_{:04d}.csv'.format(level_num))).set_index(['group','mode','cpt'])
        # Extract groups for each level
        for i_group, (group_num, sel, group_fitter) in enumerate(level):
            adp_image = os.path.join(parameterisation.out_dir, 'model/pymol/level_{}-all-modes-group_{}.png'.format(level_num, group_num))
            amp_image = os.path.join(parameterisation.out_dir, 'model/graphs/tls-model-amplitudes-level-{}-group-{}.png'.format(level_num, group_num))
            tls_mdl_str = '<br>'.join(['Mode {}:<br>T : {T11}, {T22}, {T33}, {T12}, {T13}, {T23},<br>L : {L11}, {L22}, {L33}, {L12}, {L13}, {L23},<br>S : {S11}, {S12}, {S13}, {S21}, {S22}, {S23}, {S31}, {S32}, {S33}'.format(i_mode+1, **tls_models.loc[(group_num, i_mode)].round(3)) for i_mode in xrange(parameterisation.params.fitting.tls.number_of_modes_per_group)])
            row = {
                'id'    : group_num,
                'text'  : '<br>'.join(['Number of atoms: {}'.format(sum(sel)),tls_mdl_str]),
                'adp_image' : 'data:image/png;base64,{}'.format(png2base64str(adp_image)),
                'amp_image' : 'data:image/png;base64,{}'.format(png2base64str(amp_image)),
            }
            tab['rows'].append(row)
    ##########################################################################################
    # Create tab for residual level
    ##########################################################################################
    # Create dictionary for this tab and add to tab_list
    tab = {'id'         : 'lvlres',
           'type'       : 'residual',
           'short_name' : 'Residual',
           'long_name'  : 'Final Level  (residual)',
           'description': '',
           'rows'       : [],
          }
    output_data['tabs'].append(tab)
    # Create row for each residue
    for i_rg, rg in enumerate(parameterisation.blank_master_hierarchy().select(flex.bool(parameterisation.atom_mask.tolist()),copy_atoms=True).residue_groups()):
        short_label = ShortLabeller.format(rg)
        long_label  = Labeller.format(rg)
        adp_image = os.path.join(parameterisation.out_dir, 'model/pymol/residual-residue_{}.png'.format(short_label))
        row = {
            'id'    : i_rg+1,
            'text'  : long_label,
            'adp_image' : 'data:image/png;base64,{}'.format(png2base64str(adp_image)),
        }
        tab['rows'].append(row)

    ##########################################################################################
    # Write out and format
    ##########################################################################################
    with open(out_file, 'w') as out_html:
        out_html.write(template.render(output_data).encode( "utf-8" ))
