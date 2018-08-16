import os, glob
import numpy, pandas

from scitbx.array_family import flex

from giant.structure.formatting import Labeller, ShortLabeller

from bamboo.html import png2base64src_maybe
from pandemic.html import PANDEMIC_HTML_ENV

# debugging flags
DEBUG = False

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

def wrap(string, tag='p'):
    return '<'+tag+'>'+str(string)+'</'+tag+'>'

def create_overview_tab(parameterisation):
    p = parameterisation
    f = parameterisation.fitter
    fm = p.file_manager

    chain_ids = [c.id for c in p.blank_master_hierarchy().select(flex.bool(p.atom_mask.tolist()),copy_atoms=True).chains()]

    part_f = fm.get_file('level-partitions-by-chain-template')
    prof_f = fm.get_file('png-combined-profile-template')
    resd_f = fm.get_file('png-residual-profile-template')

    tab = {'id'             : 'overview',
           'short_name'     : 'Overview',
           'long_name'      : 'Hierarchical ADP Parameterisation Overview',
           'description'    : '',
           'panels'         : [{'width':8,  'text': '', 'image':png2base64src_maybe(fm.get_file('tracking_png'), print_on_missing=DEBUG)}] + \
                              [{'width':12, 'text': ' - '.join(['<i class="fa fa-bicycle fa-fw"></i>']*25), 'class':['text-center']}] + \
                              format_summary(f.summary(), width=6) + \
                              [{'width':12, 'text': ' - '.join(['<i class="fa fa-bicycle fa-fw"></i>']*25), 'class':['text-center']}] + \
                              numpy.concatenate(zip(*[
                                  [{'width':4, 'text': 'Partition Schematic',  'image':png2base64src_maybe(part_f.format(c), print_on_missing=DEBUG)} for c in chain_ids],
                                  [{'width':4, 'text': 'TLS-level components', 'image':png2base64src_maybe(prof_f.format(c), print_on_missing=DEBUG)} for c in chain_ids],
                                  [{'width':4, 'text': 'Residual component',   'image':png2base64src_maybe(resd_f.format(c), print_on_missing=DEBUG)} for c in chain_ids],
                                                     ])).tolist() + \
                              [{'width':12, 'text': ' - '.join(['<i class="fa fa-bicycle fa-fw"></i>']*25), 'class':['text-center']}] + \
                              numpy.concatenate([format_summary(l.summary(show=False), width=6) for l in f.levels+[f.residual]]).tolist()
          }

    return tab

def create_analysis_tab(parameterisation):
    p = parameterisation
    f = parameterisation.fitter
    fm = p.file_manager

    images = []

    if p.params.analysis.calculate_r_factors:
        for f in [
                fm.get_file('all_rvalues_v_resolution'),
                fm.get_file('rgap_v_resolution'),
                fm.get_file('rfree_v_resolution'),
                fm.get_file('rwork_v_resolution'),
                ]:
            images.append({'width':6, 'image':png2base64src_maybe(f, print_on_missing=DEBUG)})

    data_table = p.tables.statistics.dropna(axis='columns', how='all')

    tab = {'id'             : 'statistics',
           'short_name'     : 'Multi-dataset Results/Analysis',
           'long_name'      : 'Model Improvement Statistics from hierarchical TLS parameterisation',
           'description'    : '',
           'plots'          : [{'div': 'variable-plots',
                                'json': data_table.T.to_json(orient='split'),
                                'default_x' : 'High Resolution Limit',
                                'default_y' : 'R-free Change (Fitted-Input)'}],
           'table'          : data_table.to_html(bold_rows=False, classes=['display nowrap'])\
                                        .replace('<th></th>','<th>Dataset</th>')\
                                        .replace('border="1" ', ''),
          }
    if images:
        tab['top_panel'] = {'heading':'R-factor Changes', 'images':images}

    return tab

def create_levels_tab(parameterisation):
    p = parameterisation
    f = parameterisation.fitter
    fm = parameterisation.file_manager

    chain_ids = [c.id for c in p.blank_master_hierarchy().select(flex.bool(p.atom_mask.tolist()),copy_atoms=True).chains()]

    tab = {'id'         : 'levels',
           'short_name' : 'Hierarchical Model Summary',
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
                    'description'   : '',
                    'panels'        : [],
                   }
    tab['tabs'].append(overview_tab)
    # Split the panels up by chain
    for c_id in chain_ids:
        # Split up the chains with divider panels
        prof_f = fm.get_file('png-combined-profile-template').format(c_id)
        resd_f = fm.get_file('png-residual-profile-template').format(c_id)
        panel = {'id'             : '<h4>Levels for Chain {}</h4>'.format(c_id),
                 'width'          : 12,
                 'show'           : True,
                 'table'          : None,
                 'contents'       : [
                     {'width':6, 'title': 'TLS-level components', 'image':png2base64src_maybe(prof_f, print_on_missing=DEBUG)},
                     {'width':6, 'title': 'Residual component',   'image':png2base64src_maybe(resd_f, print_on_missing=DEBUG)},
                                    ],
                }
        overview_tab['panels'].append(panel)
        # Add images to the overview tab for each TLS level
        for i_level, (level_num, level_lab, level) in enumerate(f):
            chain_image = fm.get_file('pml-level-chain-template').format(level_num, c_id)
            stack_image = fm.get_file('png-tls-profile-template').format(level_num, c_id)
            aniso_image = fm.get_file('png-tls-anisotropy-template').format(level_num, c_id)
            panel = {'id'             : 'Level {} of {} ({})'.format(level_num, len(f.levels),level_lab),
                     'width'          : 4 if (len(f.levels)>3) or (len(f.levels)%2==0) else 6,
                     'show'           : True,
                     'table'          : None,
                     'contents'       : [
                         {'width':12, 'image':png2base64src_maybe(chain_image, print_on_missing=DEBUG)},
                         {'width':12, 'image':png2base64src_maybe(stack_image, print_on_missing=DEBUG)},
                         {'width':12, 'image':png2base64src_maybe(aniso_image, print_on_missing=DEBUG)},
                                        ],
                    }
            overview_tab['panels'].append(panel)
        # Format residual level
        chain_image = fm.get_file('pml-residual-chain-template').format(c_id)
        stack_image = fm.get_file('png-residual-profile-template').format(c_id)
        aniso_image = fm.get_file('png-residual-anisotropy-template').format(c_id)
        panel = {'id'             : 'Final Level (residual)',
                 'width'          : 4 if (len(f.levels)>3) or (len(f.levels)%2==0) else 6,
                 'show'           : True,
                 'table'          : None,
                 'contents'       : [
                     {'width':12, 'image':png2base64src_maybe(chain_image, print_on_missing=DEBUG)},
                     {'width':12, 'image':png2base64src_maybe(stack_image, print_on_missing=DEBUG)},
                     {'width':12, 'image':png2base64src_maybe(aniso_image, print_on_missing=DEBUG)},
                                    ],
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
        # Add overview at the top of the tab
        for c_id in chain_ids:
            partn_image = fm.get_file('pml-level-partition-template').format(level_num, c_id)
            chain_image = fm.get_file('pml-level-chain-template').format(level_num, c_id)
            stack_image = fm.get_file('png-tls-profile-template').format(level_num, c_id)
            aniso_image = fm.get_file('png-tls-anisotropy-template').format(level_num, c_id)
            panel = {'id'             : 'Chain {}'.format(c_id),
                     'width'          : 12,
                     'show'           : True,
                     'table'          : None,
                     'contents'       : [
                         {'width':6, 'title':'Coloured by Group',   'image':png2base64src_maybe(partn_image, print_on_missing=DEBUG)},
                         {'width':6, 'title':'Fitted ADPs',         'image':png2base64src_maybe(chain_image, print_on_missing=DEBUG)},
                         {'width':6, 'title':'Fitted ADPs by mode', 'image':png2base64src_maybe(stack_image, print_on_missing=DEBUG)},
                         {'width':6, 'title':'Anisotropy by atom',   'image':png2base64src_maybe(aniso_image, print_on_missing=DEBUG)},
                                        ],
                    }
            level_tab['panels'].append(panel)
        # Read in the TLS models and amplitudes for this level
        tls_models     = pandas.read_csv(fm.get_file('csv-tls-mdl-template').format(level_num)).set_index(['group','model']).drop('Unnamed: 0', axis=1, errors='ignore')
        tls_amplitudes = pandas.read_csv(fm.get_file('csv-tls-amp-template').format(level_num)).set_index(['group','model','cpt']).drop('Unnamed: 0', axis=1, errors='ignore')
        # Extract groups for each level
        for i_group, (group_num, sel, group_fitter) in enumerate(level):
            # Extract TLS values for this group
            tls_vals = [tls_models.loc[(group_num, i_mode)] for i_mode in xrange(p.params.fitting.tls_models_per_tls_group)]
            # Skip if no TLS values
            if numpy.abs(tls_vals).sum() == 0.0:
                continue
            # Get images and format values
            scl_image = fm.get_file('pml-level-scaled-template').format(level_num, group_num)
            adp_image = fm.get_file('pml-level-group-template').format(level_num, group_num)
            amp_image = fm.get_file('png-tls-amp-dist-template').format(level_num, group_num)
            tls_mdl_strs = [('Mode {}:<br>' + \
                             '<samp>\n' + \
                             'T: {T11:>9.3f}, {T22:>9.3f}, {T33:>9.3f}, {T12:>9.3f}, {T13:>9.3f}, {T23:>9.3f},<br>' + \
                             'L: {L11:>9.3f}, {L22:>9.3f}, {L33:>9.3f}, {L12:>9.3f}, {L13:>9.3f}, {L23:>9.3f},<br>' + \
                             'S: {S11:>9.3f}, {S12:>9.3f}, {S13:>9.3f}, {S21:>9.3f}, {S22:>9.3f}, {S23:>9.3f}, {S31:>9.3f}, {S32:>9.3f}, {S33:>9.3f}' + \
                             '\n</samp>'
                            ).format(i_mode+1, **mode_vals.round(3)).replace(' ','&nbsp') if mode_vals.any() else 'Zero-value TLS values for mode {}'.format(i_mode+1) for i_mode, mode_vals in enumerate(tls_vals)]
            # Create panel dictionary
            panel = {'id'       : 'Group {} - {}'.format(group_num, p.levels[i_level][i_group]),
                     'width'    : 12, #max(4,12//level.n_groups()),
                     'table'    : None,
                     'contents' : [
                         {'width':12, 'text':'<br>'.join(['Number of atoms: {}'.format(sum(sel))])},
                         {'width':4,  'text':'Average size over all datasets',          'image': png2base64src_maybe(adp_image, print_on_missing=DEBUG)},
                         {'width':4,  'text':'Amplitude Distribution',                  'image': png2base64src_maybe(amp_image, print_on_missing=DEBUG)},
                         {'width':4,  'text':'Shape of disorder (arbitrary scale)',     'image': png2base64src_maybe(scl_image, print_on_missing=DEBUG)},
                                   ] + [{'width':12, 'text':s} for s in tls_mdl_strs],
                    }
            level_tab['panels'].append(panel)
        # Make  the first panel open
        if len(level_tab['panels']) > 0:
            level_tab['panels'][0]['show'] = True
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
    # Get selection for fitted atoms
    atom_sel = flex.bool(p.atom_mask.tolist())
    # Create row for each residue
    for i_chain, c in enumerate(p.blank_master_hierarchy().select(atom_sel,copy_atoms=True).chains()):
        # Panel for chain overview
        chain_image = fm.get_file('pml-residual-chain-template').format(c_id)
        stack_image = fm.get_file('png-residual-profile-template').format(c_id)
        aniso_image = fm.get_file('png-residual-anisotropy-template').format(c_id)
        panel = {'id'             : 'Residual overview for chain {}'.format(c_id),
                 'width'          : 12,
                 'show'           : True,
                 'table'          : None,
                 'contents'       : [
                     {'width':8, 'title':'Fitted ADPs',         'image':png2base64src_maybe(chain_image, print_on_missing=DEBUG)},
                     {'width':6, 'title':'Fitted ADPs profile', 'image':png2base64src_maybe(stack_image, print_on_missing=DEBUG)},
                     {'width':6, 'title':'Anisotropy by atom',   'image':png2base64src_maybe(aniso_image, print_on_missing=DEBUG)},
                                    ],
                }
        residual_tab['panels'].append(panel)
        # Panel for each residue
        panel = {'id'       : 'Residual components for chain {}'.format(c.id),
                 'width'    : 12,
                 'table'    : None,
                 'contents' : [],
                }
        residual_tab['panels'].append(panel)
        for i_rg, rg in enumerate(c.residue_groups()):
            short_label = ShortLabeller.format(rg)
            long_label  = Labeller.format(rg)
            adp_image = fm.get_file('pml-residual-group-template').format(short_label)
            panel['contents'].append({'width':4, 'text':long_label, 'image': png2base64src_maybe(adp_image, print_on_missing=DEBUG)})

    return tab

def create_settings_tab(parameterisation):

    master_phil = parameterisation.master_phil
    params = parameterisation.params
    fm = parameterisation.file_manager

    phil_str = master_phil.format(params).as_str()
    diff_str = master_phil.fetch_diff(source=master_phil.format(params)).as_str()

    tab = {'id'             : 'settings',
           'short_name'     : 'Settings',
           'long_name'      : 'Program Parameters',
           'description'    : '',
           'panels': [],
          }
    tab['panels'].append({'id'             : 'Parameters different from defaults',
                          'width'          : 12,
                          'show'           : True,
                          'contents'       : [{'width':12, 'title':'', 'text':wrap(string=diff_str, tag='pre')}],
            })
    tab['panels'].append({'id'             : 'All Parameters',
                          'width'          : 12,
                          'show'           : False,
                          'contents'       : [{'width':12, 'title':'', 'text':wrap(string=phil_str, tag='pre')}],
            })
    tab['panels'].append({'id'             : 'Parameterised penalty functions',
                          'width'          : 12,
                          'show'           : True,
                          'contents'       : [
                              {
                                  'title':'Weight of Fitting Penalty',
                                  'text' :'A function that reduces the strength of the fitting penalties as the model Uij approaches the target Uij',
                                  'image':png2base64src_maybe(fm.get_file('fitting_penalty_weights_png'), print_on_missing=DEBUG),
                                  'width':6,
                                  },
                              {
                                  'title':'Model-Input Differences (Fitting Penalty)',
                                  'text' :'Penalties for producing larger Uij that observed values',
                                  'image':png2base64src_maybe(fm.get_file('over_target_values_png'), print_on_missing=DEBUG),
                                  'width':6,
                                  },
                              {
                                  'title':'Invalid TLS Matrices (Physical Penalty)',
                                  'text' :'Penalties for producing non-physical TLS-matrices',
                                  'image':png2base64src_maybe(fm.get_file('invalid_tls_values_png'), print_on_missing=DEBUG),
                                  'width':6,
                                  },
                              {
                                  'title':'Negative Amplitudes (Physical Penalty)',
                                  'text' :'Penalties for producing negative TLS amplitudes',
                                  'image':png2base64src_maybe(fm.get_file('invalid_amplitudes_png'), print_on_missing=DEBUG),
                                  'width':6,
                                  },
                              {
                                  'title':'Invalid Uij Values (Physical Penalty)',
                                  'text' :'Penalties for producing residual Uijs with negative eigenvalues',
                                  'image':png2base64src_maybe(fm.get_file('invalid_uij_values_png'), print_on_missing=DEBUG),
                                  'width':6,
                                  },
                                            ],
            })

    return tab

def write_adp_summary(parameterisation, out_dir_tag='root'):
    """Write an overall, and a level-by-level summary of the parameterised ADPs"""

    # Get the html template
    template = PANDEMIC_HTML_ENV.get_template('adp_summary.html')
    # Create output file path
    out_file = parameterisation.file_manager.add_file(file_name='results.html', file_tag='results_html', dir_tag=out_dir_tag)

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
    # Create adp-levels tab
    ##########################################################################################
    output_data['tabs'].append(create_levels_tab(parameterisation=parameterisation))

    ##########################################################################################
    # Create statistics tab
    ##########################################################################################
    if len(parameterisation.models) > 1:
        output_data['tabs'].append(create_analysis_tab(parameterisation=parameterisation))

    ##########################################################################################
    # Create statistics tab
    ##########################################################################################
    output_data['tabs'].append(create_settings_tab(parameterisation))

    ##########################################################################################
    # Write out and format
    ##########################################################################################
    with open(out_file, 'w') as out_html:
        out_html.write(template.render(output_data).encode( "utf-8" ))
