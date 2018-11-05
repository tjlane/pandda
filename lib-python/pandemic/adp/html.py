import os, glob
import numpy, pandas

from scitbx.array_family import flex

from giant.structure.formatting import Labeller, ShortLabeller

from bamboo.html import png2base64src_maybe
from pandemic.html import PANDEMIC_HTML_ENV

# debugging flags
DEBUG = False

def format_summary(text, width=12, type='alert', colour='info'):
    paragraphs = text.split('\n\n')
    output = []
    for p in paragraphs:
        lines = p.strip('>\n ').split('\n')
        if len(lines) > 1:
            title = lines.pop(0).strip('>\n ')
        else:
            title = None
        p_dict = {'width'  : width,
                  'text'   : '<br>'.join(lines).replace('\t','&emsp;'),
                  'type'   : type,
                  'colour' : colour}
        if title is not None:
            p_dict['title'] = title
        output.append(p_dict)
    return output

def wrap(string, tag='p'):
    return '<'+tag+'>'+str(string)+'</'+tag+'>'

def create_tls_matrix_summary(tls_matrices):
    if tls_matrices.any():
        tls_mdl_str = ('<samp>' + \
                       'T: {T11:>9.3f}, {T22:>9.3f}, {T33:>9.3f}, {T12:>9.3f}, {T13:>9.3f}, {T23:>9.3f},<br>' + \
                       'L: {L11:>9.3f}, {L22:>9.3f}, {L33:>9.3f}, {L12:>9.3f}, {L13:>9.3f}, {L23:>9.3f},<br>' + \
                       'S: {S11:>9.3f}, {S12:>9.3f}, {S13:>9.3f}, {S21:>9.3f}, {S22:>9.3f}, {S23:>9.3f}, {S31:>9.3f}, {S32:>9.3f}, {S33:>9.3f}' + \
                       '</samp>').format(**tls_matrices.round(3)).replace(' ','&nbsp;')
    else:
        tls_mdl_str = "Zero-value TLS matrices"
    # Format output dictionary
    out_dict = {'text':tls_mdl_str, 'width':12, 'type':'alert', 'colour':'info'}
    return out_dict

def create_tls_decomposition_summary(tls_matrices, tolerance=1e-6):
    T = tuple(tls_matrices[['T11','T22','T33','T12','T13','T23']])
    L = tuple(tls_matrices[['L11','L22','L33','L12','L13','L23']])
    S = tuple(tls_matrices[['S11','S12','S13','S21','S22','S23','S31','S32','S33']])
    from mmtbx.tls.decompose import decompose_tls_matrices
    dcm = decompose_tls_matrices(T=T, L=L, S=S,
                                 l_and_s_in_degrees=True,
                                 tol=tolerance)
    #from IPython import embed; embed(); raise Exception()
    if dcm.is_valid():
        dcm_lines = [
                'Vibration Amplitudes: {:>9.3f}, {:>9.3f}, {:>9.3f}<br>'.format(*dcm.v_amplitudes),
                'Libration Amplitudes: {:>9.3f}, {:>9.3f}, {:>9.3f}<br>'.format(*dcm.l_amplitudes),
                'Screw Amplitudes:     {:>9.3f}, {:>9.3f}, {:>9.3f}<br>'.format(*dcm.s_amplitudes),
                '--------------------<br>',
                ('Vibration Axes:  {:>9.3f}, {:>9.3f}, {:>9.3f}<br>'+\
                 '                 {:>9.3f}, {:>9.3f}, {:>9.3f}<br>'+\
                 '                 {:>9.3f}, {:>9.3f}, {:>9.3f}<br>').format(*numpy.concatenate(dcm.v_axis_directions)),
                '--------------------<br>',
                ('Libration Axes:  {:>9.3f}, {:>9.3f}, {:>9.3f}<br>'+\
                 '                 {:>9.3f}, {:>9.3f}, {:>9.3f}<br>'+\
                 '                 {:>9.3f}, {:>9.3f}, {:>9.3f}<br>').format(*numpy.concatenate(dcm.l_axis_directions)),
                '--------------------<br>',
                ('Libration Axes   {:>9.3f}, {:>9.3f}, {:>9.3f}<br>'+\
                 'Intersections:   {:>9.3f}, {:>9.3f}, {:>9.3f}<br>'+\
                 '                 {:>9.3f}, {:>9.3f}, {:>9.3f}<br>').format(*numpy.concatenate(dcm.l_axis_intersections)),
                ]
        dcm_string = ('<samp>'+''.join(dcm_lines)+'</samp>').replace(' ','&nbsp;')
        colour = 'success'
    else:
        dcm_string = 'Unable to decompose TLS matrices: <br>&nbsp;&nbsp;&nbsp;&nbsp;{}'.format(dcm.error())
        colour = 'warning'
    # Format output dictionary
    out_dict = {'text':dcm_string, 'width':12, 'type':'alert', 'colour':colour}
    return out_dict

def create_tls_amplitude_summary(tls_amplitudes):
    tls_amp_table = tls_amplitudes.T
    # Change the names of the "Mode" columns
    new_columns = ['Mode {} - {}'.format(mode+1, cpt) for (mode, cpt) in tls_amp_table.columns]
    tls_amp_table.columns = new_columns
    # Remove columns names and add index name
    tls_amp_table.columns.name = 'Dataset'
    tls_amp_table = tls_amp_table.to_html(bold_rows=False, classes=['table table-condensed table-hover datatable nowrap text-center'], border=0)
    tls_amp_table = tls_amp_table.replace('<th>', '<th class="text-center">')
    return {'table': tls_amp_table}

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
           'panels'         : [
               {
                   'title'  : 'Parameterisation Summary',
                   'body'   : [{'width':8,  'text': '', 'image':png2base64src_maybe(fm.get_file('tracking_png'), print_on_missing=DEBUG)}] + \
                              [{'width':12, 'text': ' - '.join(['<i class="fa fa-bicycle fa-fw"></i>']*25), 'classes':['text-center']}] + \
                              format_summary(f.summary(), width=6),
                   },
               {
                   'title':'Chain-by-Chain Disorder Summaries',
                   'body': numpy.concatenate(zip(*[
                                   [{'width':12,'title': 'Chain {}'.format(c)} for c in chain_ids],
                                   [{'width':4, 'text': 'Partition Schematic',  'image':png2base64src_maybe(part_f.format(c), print_on_missing=DEBUG)} for c in chain_ids],
                                   [{'width':4, 'text': 'TLS-level components', 'image':png2base64src_maybe(prof_f.format(c), print_on_missing=DEBUG)} for c in chain_ids],
                                   [{'width':4, 'text': 'Residual component',   'image':png2base64src_maybe(resd_f.format(c), print_on_missing=DEBUG)} for c in chain_ids],
                                            ])).tolist(),
                   },
               {
                   'title' : 'Level-by-Level Model Parameters Summary',
                   'body'  : numpy.concatenate([format_summary(l.summary(show=False), width=6) for l in f.levels+[f.residual]]).tolist(),
                   },
               ],
          }

    return tab

def create_levels_tab(parameterisation):
    p = parameterisation
    f = parameterisation.fitter
    fm = parameterisation.file_manager

    tls_tolerance = p.params.fitting.precision.tls_tolerance

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
        panel = {'title' : '<h4>Levels for Chain {}</h4>'.format(c_id),
                 'body'  : [
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
            panel = {'title' : 'Level {} of {} ({})'.format(level_num, len(f.levels),level_lab),
                     'width' : 4 if (len(f.levels)>3) or (len(f.levels)%2==0) else 6,
                     'body'  : [
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
        panel = {'title' : 'Final Level (residual)',
                 'width' : 4 if (len(f.levels)>3) or (len(f.levels)%2==0) else 6,
                 'body'  : [
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
            panel = {'title' : 'Chain {}'.format(c_id),
                     'body'  : [
                         {'width':6, 'title':'Coloured by Group',   'image':png2base64src_maybe(partn_image, print_on_missing=DEBUG)},
                         {'width':6, 'title':'Fitted ADPs',         'image':png2base64src_maybe(chain_image, print_on_missing=DEBUG)},
                         {'width':6, 'title':'Fitted ADPs by mode', 'image':png2base64src_maybe(stack_image, print_on_missing=DEBUG)},
                         {'width':6, 'title':'Anisotropy by atom',  'image':png2base64src_maybe(aniso_image, print_on_missing=DEBUG)},
                                        ],
                    }
            level_tab['panels'].append(panel)
        # Read in the TLS models and amplitudes for this level
        tls_models     = pandas.read_csv(fm.get_file('csv-tls-mdl-template').format(level_num)).set_index(['group','model']).drop('Unnamed: 0', axis=1, errors='ignore')
        tls_amplitudes = pandas.read_csv(fm.get_file('csv-tls-amp-template').format(level_num)).set_index(['group','model','cpt']).drop('Unnamed: 0', axis=1, errors='ignore')
        # Extract groups for each level
        for i_group, (group_num, sel, group_fitter) in enumerate(level):
            # Extract TLS values for this group
            tls_mats = tls_models.loc[group_num]
            tls_amps = tls_amplitudes.loc[group_num]
            # Get TLS matrix summaries
            tls_mat_dicts = [create_tls_matrix_summary(tls_matrices=mats) for idx, mats in tls_mats.iterrows()]
            for i, d in enumerate(tls_mat_dicts):
                d.update({'title':'Mode {}:'.format(i+1)})
            tls_dcm_dicts = [create_tls_decomposition_summary(tls_matrices=mats, tolerance=tls_tolerance) for idx, mats in tls_mats.iterrows()]
            dcm_width = 12 if (p.params.fitting.tls_models_per_tls_group==1) else 6
            for i, d in enumerate(tls_dcm_dicts):
                d.update({'title':'TLS Decomposition of Mode {}:'.format(i+1), 'width':dcm_width})
            # Get TLS amplitude summaries
            tls_amp_dict = create_tls_amplitude_summary(tls_amplitudes=tls_amps)
            tls_amp_dict.update({'width':12, 'type':'alert', 'colour':'info', 'title':'All Amplitudes'})
            # Get images and format values
            scl_image = fm.get_file('pml-level-scaled-template').format(level_num, group_num)
            adp_image = fm.get_file('pml-level-group-template').format(level_num, group_num)
            amp_image = fm.get_file('png-tls-amp-dist-template').format(level_num, group_num)
            image_dicts = []
            for text, image in [
                    ('Average size over all datasets',      adp_image),
                    ('Amplitude Distribution',              amp_image),
                    ('Shape of disorder (arbitrary scale)', scl_image),
                                ]:
                if os.path.exists(image):
                    image_dicts.append({'width':4, 'text':text, 'image':png2base64src_maybe(image, print_on_missing=DEBUG)})
            # Create panel dictionary
            panel = {'title' : 'Group {} - {}'.format(group_num, p.levels[i_level][i_group]),
                     'width' : 12, #max(4,12//level.n_groups()),
                     'show'  : (i_group==0),
                     'body'  : [
                         {'width':12, 'text':'<br>'.join(['Number of atoms: {}'.format(sum(sel))])}
                         ] + tls_mat_dicts + image_dicts + tls_dcm_dicts + [tls_amp_dict]
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
        panel = {'title' : 'Residual overview for chain {}'.format(c_id),
                 'body'  : [
                     {'width':8, 'title':'Fitted ADPs',         'image':png2base64src_maybe(chain_image, print_on_missing=DEBUG)},
                     {'width':6, 'title':'Fitted ADPs profile', 'image':png2base64src_maybe(stack_image, print_on_missing=DEBUG)},
                     {'width':6, 'title':'Anisotropy by atom',   'image':png2base64src_maybe(aniso_image, print_on_missing=DEBUG)},
                                    ],
                }
        residual_tab['panels'].append(panel)
        # Panel for each residue
        panel = {'title' : 'Residual components for chain {}'.format(c.id),
                 'show'  : False,
                 'body'  : [],
                }
        residual_tab['panels'].append(panel)
        for i_rg, rg in enumerate(c.residue_groups()):
            short_label = ShortLabeller.format(rg)
            long_label  = Labeller.format(rg)
            panel['body'].append({'width':4, 'text':long_label})
            adp_image = fm.get_file('pml-residual-group-template').format(short_label)
            if os.path.exists(adp_image):
                panel['body'][-1]['image'] = png2base64src_maybe(adp_image, print_on_missing=DEBUG)

    return tab

def create_analysis_tab(parameterisation):
    p = parameterisation
    f = parameterisation.fitter
    fm = p.file_manager

    tab = {'id'             : 'statistics',
           'short_name'     : 'Multi-dataset Results/Analysis',
           'long_name'      : 'Model Improvement Statistics from hierarchical TLS parameterisation',
           'description'    : '',
           'panels'         : [],
           }

    # Add panel containing R-factor plots images
    if p.params.analysis.calculate_r_factors:
        images = []
        for f in [fm.get_file('all_rvalues_v_resolution'),
                  fm.get_file('rgap_v_resolution'),
                  fm.get_file('rfree_v_resolution'),
                  fm.get_file('rwork_v_resolution'),
                  ]:
            images.append({'width':6, 'image':png2base64src_maybe(f, print_on_missing=DEBUG)})
        # Append to panels
        tab['panels'].append({
            'title' :'R-factor Changes',
            'body'  : images
            })

    # Extract data table for plot and table
    data_table = p.tables.statistics.dropna(axis='columns', how='all')
    # Panel for the interactive plots
    json_plot = {'id'        : 'variable-plots',
                 'json'      : data_table.T.to_json(orient='split'),
                 'default_x' : 'High Resolution Limit',
                 'default_y' : 'R-free Change (Fitted-Input)'}
    tab['panels'].append({
        'title' : 'Interactive Summary Graphs',
        'body'  : [{'divs': [json_plot]}],
        })
    # Panel for the table of data
    tab['panels'].append({
        'title' : 'Inidividual Dataset Statistics',
        'body'  : [{'text': 'Data from output CSV (dataset_scores.csv)'},
                   {'table': data_table.to_html(bold_rows=False, classes=['table table-striped table-hover datatable nowrap'])\
                           .replace('<th></th>','<th>Dataset</th>')\
                           .replace('border="1" ', '')},
                   ],
        })

    return tab, [json_plot]

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
    tab['panels'].append({'title'          : 'Parameters different from defaults',
                          'width'          : 12,
                          'body'           : [{'width':12, 'title':'', 'text':wrap(string=diff_str, tag='pre')}],
            })
    tab['panels'].append({'title'          : 'All Parameters',
                          'width'          : 12,
                          'show'           : False,
                          'body'           : [{'width':12, 'title':'', 'text':wrap(string=phil_str, tag='pre')}],
            })
    tab['panels'].append({'title'          : 'Parameterised penalty functions',
                          'width'          : 12,
                          'body'           : [
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
        tab, json_plots = create_analysis_tab(parameterisation=parameterisation)
        output_data['tabs'].append(tab)
        output_data.setdefault('json_plots',[]).extend(json_plots)

    ##########################################################################################
    # Create statistics tab
    ##########################################################################################
    output_data['tabs'].append(create_settings_tab(parameterisation))

    ##########################################################################################
    # Write out and format
    ##########################################################################################
    with open(out_file, 'w') as out_html:
        out_html.write(template.render(output_data).encode( "utf-8" ))
