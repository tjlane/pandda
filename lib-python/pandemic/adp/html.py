import os, glob
import numpy, pandas

from scitbx.array_family import flex

from giant.structure.formatting import Labeller, ShortLabeller

from bamboo.html import png2base64src
from pandemic.html import PANDEMIC_HTML_ENV

def format_summary(text, width=12, cls_tuple=['alert','info']):
    paragraphs = text.split('\n\n')
    output = []
    for p in paragraphs:
        lines = p.strip('>\n ').split('\n')
        if len(lines) > 1:
            title = lines.pop(0).strip('>\n ')
        else:
            title = None
        p_dict = {'width' : width,
                  'text'  : '<br>'.join(lines).replace('\t','&emsp;'),
                  'class' : cls_tuple}
        if title is not None:
            p_dict['title'] = title
        output.append(p_dict)
    return output

def create_overview_tab(parameterisation):
    p = parameterisation
    f = parameterisation.fitter

    chain_ids = [c.id for c in p.blank_master_hierarchy().select(flex.bool(p.atom_mask.tolist()),copy_atoms=True).chains()]

    tab = {'id'             : 'overview',
           'short_name'     : 'Overview',
           'long_name'      : 'Hierarchical ADP Parameterisation Overview',
           'description'    : '',
           'contents'       : format_summary(f.summary(), width=6) + \
                              numpy.concatenate([format_summary(l.summary(show=False), width=3) for l in f.levels]).tolist() + \
                              #[format_summary(f.residual.summary(), width=12)] + \
                              [{'width':12, 'text': '---'}] + \
                              # [{'width': 12,
                              #  'title': 'Hierarchy Summary:',
                              #  'text':  f.html_summary(),#'{} levels (+residual level):{}'.format(len(p.levels), '<br>&nbsp;'.join(['']+['Level {}: {}'.format(il+1, l) for il,l in enumerate(p.level_labels)])),
                              #  'class': ['alert','info'],
                              # }] +\
                              numpy.concatenate(zip(*[[{'width':4, 'text': 'Partition Schematic',  'image':png2base64src(os.path.join(p.out_dir, 'model/level-partitioning-chain-{}.png'.format(c)))    } for c in chain_ids],
                                                      [{'width':4, 'text': 'TLS-level components', 'image':png2base64src(os.path.join(p.out_dir, 'model/all-stacked-chain_{}.png'.format(c)))           } for c in chain_ids],
                                                      [{'width':4, 'text': 'Residual component',   'image':png2base64src(os.path.join(p.out_dir, 'model/all-residual-chain_{}.png'.format(c)))          } for c in chain_ids]
                                                     ])).tolist(),
          }

    return tab

def create_analysis_tab(parameterisation):
    p = parameterisation
    f = parameterisation.fitter

    tab = {'id'             : 'statistics',
           'short_name'     : 'Results/Analysis',
           'long_name'      : 'Model Improvement Statistics from hierarchical TLS parameterisation',
           'description'    : 'things',
           'images'         : [{'width':6, 'path':png2base64src(os.path.join(p.out_dir, 'graphs/r-values-change-by-resolution.png'))},
                               {'width':6, 'path':png2base64src(os.path.join(p.out_dir, 'graphs/r-gap-by-resolution.png'))},
                               {'width':6, 'path':png2base64src(os.path.join(p.out_dir, 'graphs/r-free-by-resolution.png'))},
                               {'width':6, 'path':png2base64src(os.path.join(p.out_dir, 'graphs/r-work-by-resolution.png'))}],
           'plots'          : [{'div': 'variable-plots',
                                'json': p.table.T.to_json(orient='split'),
                                'default_x' : 'High Resolution Limit',
                                'default_y' : 'R-free Change (Fitted-Input)'}],
           'table'          : p.table.to_html(bold_rows=False, classes=['display nowrap'])\
                                        .replace('<th></th>','<th>Dataset</th>')\
                                        .replace('border="1" ', ''),
          }

    return tab

def create_levels_tab(parameterisation):
    p = parameterisation
    f = parameterisation.fitter

    tab = {'id'         : 'levels',
           'short_name' : 'ADP Summary',
           'long_name'  : 'Level-by-level TLS parameterisation',
           'description': 'Parameteristaion composed of {} levels.'.format(len(f.levels)),
           'tabs'       : [],
          }
    # -------------------------------->
    # Create overview sub-tab
    # -------------------------------->
    overview_tab = {'id'            : tab['id']+'overview',
                    'active'        : True,
                    'short_name'    : 'Overview',
                    'long_name'     : 'Overview of the parameterised hierarchical ADP model',
                    'description'   : 'things',
                    'panels'        : [],
                   }
    tab['tabs'].append(overview_tab)
    # Add images to the overview tab for each TLS level
    for i_level, (level_num, level_lab, level) in enumerate(f):
        chain_image = os.path.join(p.out_dir, 'model/pymol/level_{}-all-modes-chain_A.png'.format(level_num))
        stack_image = os.path.join(p.out_dir, 'model/level_{}-all-modes-TLS-chain_A.png'.format(level_num))
        panel = {'id'             : 'Level {} of {} ({})'.format(level_num, len(f.levels),level_lab),
                 'text'           : '{} atoms.'.format('X'),
                 'width'          : 4,
                 'table'          : None,
                 'images'         : [{'width':12, 'path':png2base64src(chain_image) if os.path.exists(chain_image) else ''},
                                     {'width':12, 'path':png2base64src(stack_image) if os.path.exists(stack_image) else ''}],
                }
        overview_tab['panels'].append(panel)
    # Format residual level
    chain_image = os.path.join(p.out_dir, 'model/pymol/residual-chain_A.png')
    stack_image = os.path.join(p.out_dir, 'model/all-residual-chain_A.png')
    panel = {'id'             : 'Final Level (residual)',
             'text'           : '{} atoms.'.format('X'),
             'width'          : 4,
             'table'          : None,
             'images'         : [{'width':12, 'path':png2base64src(chain_image) if os.path.exists(chain_image) else ''},
                                 {'width':12, 'path':png2base64src(stack_image) if os.path.exists(stack_image) else ''}],
              }
    overview_tab['panels'].append(panel)
    # -------------------------------->
    # Create tab for each level
    # -------------------------------->
    for i_level, (level_num, level_lab, level) in enumerate(f):
        # Create dictionary for this tab and add to tab_list
        level_tab = {'id'         : tab['id']+'lvl{}'.format(level_num),
                     'short_name' : 'Level {}'.format(level_num),
                     'long_name'  : 'Level {} ({})'.format(level_num, level_lab),
                     'description': 'Level {} of {}. '.format(level_num, len(f.levels))+\
                                    'Composed of {} groups'.format(level.n_groups()),
                     'panels'       : [],
              }
        tab['tabs'].append(level_tab)
        # Read in the TLS models and amplitudes for this level
        tls_models     = pandas.read_csv(os.path.join(p.out_dir, 'model/csvs/tls_models_level_{:04d}.csv'.format(level_num))).set_index(['group','mode'])
        tls_amplitudes = pandas.read_csv(os.path.join(p.out_dir, 'model/csvs/tls_amplitudes_level_{:04d}.csv'.format(level_num))).set_index(['group','mode','cpt'])
        # Extract groups for each level
        for i_group, (group_num, sel, group_fitter) in enumerate(level):
            adp_image = os.path.join(p.out_dir, 'model/pymol/level_{}-all-modes-group_{}.png'.format(level_num, group_num))
            amp_image = os.path.join(p.out_dir, 'model/graphs/tls-model-amplitudes-level-{}-group-{}.png'.format(level_num, group_num))
            tls_mdl_str = '<br>'.join(['Mode {}:<br>T : {T11}, {T22}, {T33}, {T12}, {T13}, {T23},<br>L : {L11}, {L22}, {L33}, {L12}, {L13}, {L23},<br>S : {S11}, {S12}, {S13}, {S21}, {S22}, {S23}, {S31}, {S32}, {S33}'.format(i_mode+1, **tls_models.loc[(group_num, i_mode)].round(3)) for i_mode in xrange(p.params.fitting.number_of_modes_per_group)])
            panel = {'id'    : 'Group {} - {}'.format(group_num, p.levels[i_level][i_group]),
                     'text'  : '<br>'.join(['Number of atoms: {}'.format(sum(sel)),tls_mdl_str]),
                     'width' : max(4,12//level.n_groups()),
                     'table' : None,
                     'images': [{'width':12, 'path': png2base64src(adp_image)},
                                {'width':12, 'path': png2base64src(amp_image)}],
                    }
            level_tab['panels'].append(panel)
    # -------------------------------->
    # Create tab for residual level
    # -------------------------------->
    residual_tab = {'id'         : 'lvlres',
                    'short_name' : 'Residual',
                    'long_name'  : 'Final Level  (residual)',
                    'description': '',
                    'panels'     : [],
                   }
    tab['tabs'].append(residual_tab)
    # Create row for each residue
    for i_rg, rg in enumerate(p.blank_master_hierarchy().select(flex.bool(p.atom_mask.tolist()),copy_atoms=True).residue_groups()):
        short_label = ShortLabeller.format(rg)
        long_label  = Labeller.format(rg)
        adp_image = os.path.join(p.out_dir, 'model/pymol/residual-residue_{}.png'.format(short_label))
        panel = {'id'    : long_label,
                 'width' : 4,
                 'text'  : '',
                 'table' : None,
                 'images': [{'width':12, 'path': png2base64src(adp_image)}],
                }
        residual_tab['panels'].append(panel)

    return tab

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
    output_data['header'] = 'PanDEMIC ADP Output Summary'
    output_data['title'] = 'PanDEMIC ADP Parameterisation'
    output_data['introduction'] = 'Hierarchical ADP parameterisation summary from pandemic.adp'
    # ===========================================================>
    # Construct the tabs
    output_data['tabs'] = []

    ##########################################################################################
    # Create overview tab
    ##########################################################################################
    output_data['tabs'].append(create_overview_tab(parameterisation=parameterisation))
    output_data['tabs'][-1]['active'] = True
    ##########################################################################################
    # Create statistics tab
    ##########################################################################################
    output_data['tabs'].append(create_analysis_tab(parameterisation=parameterisation))
    ##########################################################################################
    # Create adp-levels tab
    ##########################################################################################
    output_data['tabs'].append(create_levels_tab(parameterisation=parameterisation))
    ##########################################################################################
    # Write out and format
    ##########################################################################################
    with open(out_file, 'w') as out_html:
        out_html.write(template.render(output_data).encode( "utf-8" ))
