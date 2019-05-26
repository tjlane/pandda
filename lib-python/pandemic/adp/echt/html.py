from libtbx import adopt_init_args

import math
import numpy, pandas

# debugging flags
DEBUG = False


from pandemic.adp.html.parameters import ParameterHtmlSummary
class EchtParameterHtmlSummary(ParameterHtmlSummary):


    def main_summary(self):

        pf = self.parameter_files

        panels = []
        if pf.get('barrier_penalty'):
            panels.append({
                'type'      : 'panel',
                'title'     : 'Simplex optimisation functions',
                'width'     : 12,
                'contents'  : [
                    {
                        'width' : 6,
                        'title' : 'Model-Input Differences (Fitting Penalty)',
                        'text'  : 'Penalties for producing larger Uij that observed values (used in the first cycle only)',
                        'image' : self.image(pf.get('barrier_penalty')),
                        },
                    ],
                })

        if pf.get('level amplitudes weights'):
            wgt_dict = pf.get('level amplitudes weights')
            panels.append({
                'type'      : 'panel',
                'title'     : 'Level amplitude optimisation weights',
                'width'     : 12,
                'contents'  : [{
                        'width' : 6,
                        'title' : variable,
                        'image' : self.image(image_path),
                        } for variable, image_path in wgt_dict.iteritems()],
                })

        self.main_output.setdefault('contents', []).extend(panels)
        return [self.main_output]


from pandemic.adp.html import HtmlSummary
class EchtModelHtmlSummary(HtmlSummary):


    DEG2RAD = math.pi/180.0
    DEG2RADSQ = DEG2RAD*DEG2RAD
    RAD2DEG = 1.0 / DEG2RAD
    RAD2DEGSQ = 1.0 / DEG2RADSQ

    def __init__(self,
        hierarchy_files,
        model_files,
        model_object,
        isotropic_mask,
        parameters,
        ):
        adopt_init_args(self, locals())

    def short_summary(self):

        hf = self.hierarchy_files
        mf = self.model_files

        img1 = mf.get('all_levels_uijs_profiles_png', {})
        img2 = mf.get('all_levels_uijs_anisotropy_png', {})
        img3 = hf.get('level_partitions_png', {})

        chain_ids = sorted(set(img1.keys()+img2.keys()+img3.keys()))

        output = {
            'type'     : 'panel',
            'title'    : 'Chain-by-Chain Disorder Summaries',
            'contents' : [],
            }

        for c_id in chain_ids:
            output.setdefault('contents', []).extend([
                {
                    'width' : 12,
                    'title' : 'Chain {}'.format(c_id),
                    },
                {
                    'width' : 4,
                    'text'  : 'Level-by-level ADP Profile',
                    'image' : self.image(img1.get(c_id)),
                    },
                {
                    'width' : 4,
                    'text'  : 'Total model anisotropy',
                    'image' : self.image(img2.get(c_id)),
                    },
                {
                    'width' : 4,
                    'text'  : 'Partition Schematic',
                    'image' : self.image(img3.get(c_id)),
                    },
                ])

        return [output]

    def main_summary(self):

        tab = {'id'        : 'levels',
               'alt_title' : 'Hierarchical Model',
               'title'     : 'ECHT level-by-level TLS parameterisation',
               'contents': [],
              }

        tab['contents'] += self.model_hierarchy_summary_string_list(
            model_object = self.model_object,
            isotropic_mask = self.isotropic_mask,
            )
        tab['contents'] += self.model_parameter_summary_string_list(
            model_object = self.model_object,
            isotropic_mask = self.isotropic_mask,
            )

        overview_tab = self.overview_tab(parent_tab_id=tab['id'])
        level_tabs = self.create_tls_level_tabs(parent_tab_id=tab['id'])

        tab_set = {
            'type' : 'tabs',
            'contents' : [overview_tab] + level_tabs
        }

        tab.setdefault('contents', []).append(tab_set)

        return [tab]

    def model_hierarchy_summary_string_list(self,
        model_object,
        isotropic_mask,
        ):

        mo = model_object

        n_tls_groups = numpy.sum([len(gs) for gs in mo.tls_objects])

        atoms_per_group_per_level = [[g.n_atoms for g in gs] for gs in mo.tls_objects]
        group_sizes_per_level = ['{} - {} atoms'.format(min(n), max(n)) for n in atoms_per_group_per_level]

        if (isotropic_mask is not None):
            n_iso = isotropic_mask.selection.iselection().size()
            n_ani = mo.n_atoms - n_iso
        else:
            n_iso = 0
            n_ani = mo.n_atoms

        adp_level_opt = numpy.array(mo.adp_values).any()

        s = ''
        s += '\n> Hierarchical Model Summary:'
        s += '\nNumber of TLS modes per group: {}'.format(mo.n_modes)
        s += '\nModel has {} TLS Levels:'.format(mo.n_tls_levels)
        for i_l, l in enumerate(mo.tls_level_names):
            s += '\n\tLevel {} ({}): {} groups (sizes {})'.format(
                i_l+1,
                l,
                len(mo.tls_objects[i_l]),
                group_sizes_per_level[i_l],
                )
        s += '\nTotal TLS Groups: {}'.format(n_tls_groups)
        if adp_level_opt:
            s += '\nSize of ADP level ({}): {} atoms'.format(mo.adp_level_name, mo.n_atoms)
            s += '\nAnisotropic - Isotropic atoms: {} - {}'.format(n_ani, n_iso)
        else:
            s += '\n\tNo atomic level was optimised'.format(mo.adp_level_name)

        return self.format_summary(s, width=6)

    def model_parameter_summary_string_list(self,
        model_object,
        isotropic_mask,
        ):

        mo = model_object

        n_tls_params_per_group = mo.n_modes*(20+mo.n_datasets-1)
        n_tls_groups_per_level = [len(gs) for gs in mo.tls_objects]
        n_tls_params_per_level = [n_tls_params_per_group*n for n in n_tls_groups_per_level]

        adp_level_opt = numpy.array(mo.adp_values).any()

        if (not adp_level_opt):
            n_iso = 0
            n_ani = 0
        elif (isotropic_mask is not None):
            n_iso = isotropic_mask.selection.iselection().size()
            n_ani = isotropic_mask.selection.size() - n_iso
        else:
            n_iso = 0
            n_ani = mo.n_atoms

        n_adp_params = 6 * n_ani + n_iso
        n_tls_params = sum(n_tls_params_per_level)
        n_params = n_adp_params + n_tls_params

        s = ''
        s += '\n> Parameterisation Summary:'
        s += '\nParameters per group: {} (20 from matrices + {} from amplitudes)'.format(n_tls_params_per_group, mo.n_datasets-1)
        s += '\nParameters per TLS level:'
        for i_l, l in enumerate(mo.tls_level_names):
            s += '\n\tLevel {} ({}): {}'.format(
                i_l+1,
                l,
                n_tls_params_per_level[i_l],
                )
        s += '\nTotal TLS parameters: {}'.format(n_tls_params)
        if adp_level_opt:
            s += '\nTotal ADP ({} level) parameters: {} (6*{} + 1*{})'.format(mo.adp_level_name, n_adp_params, n_ani, n_iso)
        s += '\nTotal parameters: {} ({:.2f} per atom per dataset)'.format(n_params, n_params/(mo.n_atoms * mo.n_datasets))

        return self.format_summary(s, width=6)

    def overview_tab(self, parent_tab_id):
        """Make overview tab for the hierarchical model"""

        hf = self.hierarchy_files   # definitions
        mf = self.model_files       # fitted

        img1 = mf.get('all_levels_uijs_profiles_png', {})
        img2 = mf.get('all_levels_uijs_anisotropy_png', {})

        img3 = mf.get('level_uijs_pymol_by_chain_png', {})
        img4 = mf.get('level_uijs_profiles_png', {})
        img5 = mf.get('level_uijs_anisotropy_png', {})

        chain_ids = sorted(img1.keys())

        # -------------------------------->
        # Create overview sub-tab
        # -------------------------------->
        overview_tab = {'id'        : parent_tab_id+'overview',
                        'active'    : True,
                        'alt_title' : 'Overview',
                        'title'     : 'Overview of the parameterised hierarchical ADP model',
                        'contents'  : [],
                       }

        tabs = {
            'type' : 'tabs',
            'contents' : [],
            }
        overview_tab.setdefault('contents', []).append(tabs)

        # Split the panels up by chain
        for c_id in chain_ids:

            # Split up the chains with divider panels
            tab = {
                'alt_title' : 'Chain {}'.format(c_id),
                'title' : 'Levels for Chain {}'.format(c_id),
                'contents'  : [
                    {
                        'width' : 6,
                        'title' : 'Level-by-level Disorder Profile',
                        'image' : self.image(img1.get(c_id)),
                        },
                    {
                        'width' : 6,
                        'title' : 'Anisotropy of Total Disorder',
                        'image' : self.image(img2.get(c_id)),
                        },
                    ],
                }
            tabs.setdefault('contents', []).append(tab)

            n_levels = self.model_object.n_levels
            for i_l, l in enumerate(self.model_object.all_level_names):

                d = {
                        'title' : 'Level {} of {} ({})'.format(i_l+1, n_levels, l),
                        'width' : 12, #4 if (n_levels>3) or (n_levels%2==0) else 6
                        'type'  : 'alert',
                        'colour': 'info',
                        'contents'  : [
                            {
                                'width' : 4,
                                'image' : self.image(img3.get(l, {}).get(c_id)),
                                },
                            {
                                'width' : 4,
                                'image' : self.image(img4.get(l, {}).get(c_id)),
                                },
                            {
                                'width' : 4,
                                'image' : self.image(img5.get(l, {}).get(c_id)),
                                },
                            ],
                        }
                tab.setdefault('contents', []).append(d)

        tabs['contents'][0].setdefault('active', True)

        return overview_tab

    def create_tls_level_tabs(self, parent_tab_id):

        hf = self.hierarchy_files
        mf = self.model_files

        mo = self.model_object

        img1 = hf.get('level_atoms_pymol_png', {})
        img2 = mf.get('level_uijs_pymol_by_chain_png', {})
        img3 = mf.get('level_uijs_profiles_png', {})
        img4 = mf.get('level_uijs_anisotropy_png', {})

        chain_ids = sorted(img1.values()[0].keys())

        # -------------------------------->
        # Create tab for each level
        # -------------------------------->
        tabs = []
        for i_l, l in enumerate(mo.tls_level_names):

            # Create dictionary for this tab and add to tab_list
            level_tab = {'id'         : parent_tab_id+'lvl{}'.format(i_l),
                         'alt_title' : 'Level {}'.format(i_l+1),
                         'title'  : 'Level {} ({})'.format(i_l+1, l),
                         'contents'   : [
                            {
                                'text' : 'Level {} of {}. '.format(i_l+1, mo.n_levels)+\
                                         'Composed of {} groups'.format(len(mo.tls_objects[i_l])),
                                },
                         ],
                  }
            tabs.append(level_tab)

            ########################################################

            # Add overview at the top of the tab
            tab_set = {
                'type'  : 'tabs',
                'title' : "Chain-by-Chain Summaries",
                'contents'  : [],
            }
            level_tab.setdefault('contents',[]).append(tab_set)
            for c_id in chain_ids:

                chain_tab = {
                    'id'        : level_tab['id']+'chain{}'.format(self.counter.next()),
                    'alt_title' : 'Chain {}'.format(c_id),
                    'title'     : 'Summary for Chain {}'.format(c_id),
                    'contents'  : [
                        {
                            'width' : 6,
                            'title' : 'Image of chain, coloured by group',
                            'image' : self.image(img1.get(l,{}).get(c_id)),
                            },
                        {
                            'width' : 6,
                            'title' : 'Level ADPs' + ' (average size over all datasets)'*(mo.n_datasets>1),
                            'image' : self.image(img2.get(l,{}).get(c_id)),
                            'body'  : [{'text' :'<strong>Please note</strong>: highly anisotropic atoms may not be visible in the above image.'}],
                            },
                        {
                            'width' : 6,
                            'title' : 'Level ADPs by mode' + ' (average size over all datasets)'*(mo.n_datasets>1),
                            'image' : self.image(img3.get(l,{}).get(c_id)),
                            },
                        {
                            'width' : 6,
                            'title' : 'ADP Anisotropy' + ' (from average over all datasets)'*(mo.n_datasets>1),
                            'image' : self.image(img4.get(l,{}).get(c_id)),
                            },
                        ],
                    }
                tab_set.setdefault('contents',[]).append(chain_tab)
            tab_set['contents'][0].setdefault('active', True)

            ########################################################

            tab_set = {
                'type'  : 'tabs',
                'title' : "Group-by-Group Summaries",
                'classes' : ['nav-hover'],
                'contents'  : [],
            }
            level_tab.setdefault('contents',[]).append(tab_set)

            # Read in the TLS modes and amplitudes for this level
            tls_matrices   = pandas.read_csv(mf['level_tls_matrices_csv'][l]).set_index(['group','mode']).drop('Unnamed: 0', axis=1, errors='ignore')
            tls_amplitudes = pandas.read_csv(mf['level_tls_amplitudes_csv'][l]).set_index(['group','mode']).drop('Unnamed: 0', axis=1, errors='ignore')

            # Extract groups for each level
            for i_g, tls_obj in enumerate(mo.tls_objects[i_l]):
                n_g = i_g + 1

                # Extract TLS values for this group
                tls_mats = tls_matrices.loc[n_g].drop(columns=['label'])
                tls_amps = tls_amplitudes.loc[n_g].drop(columns=['label'])
                tls_amps_av = tls_amps.mean(axis=1).values # average over datasets
                tls_coms = pandas.DataFrame(index=tls_amps.T.index, data=numpy.array(tls_obj.origins))

                # Get TLS matrix summaries for core model
                dicts = [{'title':'TLS matrices'}]
                for i_mode, (idx, mats) in enumerate(tls_mats.iterrows()):
                    d = self.create_tls_matrix_summary(tls_matrices=mats)
                    d.update({
                        'width' : 12,
                        'title' : 'Mode {} (normalised to arbitrary scale)'.format(idx),
                        })
                    dicts.append(d)
                tls_mat_dicts_core = dicts

                # Get TLS matrix summaries for each dataset
                tls_mat_dicts = []
                for i_dst, (dst_lab, dst_amps) in enumerate(tls_amps.T.iterrows()):

                    dicts = []
                    d = {
                            'type' : 'alert',
                            'colour' : 'info',
                            'title' : 'TLS amplitudes',
                            'text' : '<br>'.join(['<samp>Mode {}: <strong>{:>7.5f}</strong></samp>'.format(i+1, v) for i,v in enumerate(dst_amps.values)]),
                        }
                    dicts.append(d)

                    for i_mode, (idx, mats) in enumerate(tls_mats.iterrows()):
                        # TLS matrix summary
                        d = self.create_tls_matrix_summary(tls_matrices=dst_amps.values[i_mode]*mats)
                        d.update({
                            'width' : 12,
                            'title' : 'Mode {} (scaled for this dataset):'.format(i_mode+1),
                            })
                        dicts.append(d)
                        # TLS decomposition summary
                        d = self.create_tls_decomposition_summary(
                            tls_matrices = mats,
                            tolerance = self.parameters.optimisation.tolerances.tls_matrix_tolerance,
                            amplitude = dst_amps.values[i_mode],
                            )
                        d.update({
                            'width' : 12,
                            'title' : 'TLS decomposition of mode {}:'.format(i_mode+1),
                            })
                        dicts.append(d)

                    tls_mat_dicts += [
                            {
                                'title' : 'TLS matrices for dataset {}'.format(dst_lab),
                                'classes' : ['text-center'],
                                'text' :
                                    '<i class="fa fa-angle-double-left fa-2x" style="vertical-align: middle;"></i>' + \
                                    ' <strong>scroll for other datasets</strong> ' + \
                                    '<i class="fa fa-angle-double-right fa-2x" style="vertical-align: middle;"></i>',
                                'width' : 12,
                                'contents' : dicts,
                                },
                            ]

                # Create scrolling div to contain summaries
                tls_mat_dicts_dsts = {
                    'contents' : [
                        {
                            'type'  : 'panel',
                            'title' : 'TLS matrices for each dataset (click to expand/collapse)',
                            'show'  : True,
                            'contents':[
                                {
                                    'type' : 'scroll',
                                    'contents' : tls_mat_dicts,
                                    },
                                ],
                            },
                        ],
                    }

                # Get TLS amplitude summaries
                tls_amp_dict = self.create_tls_amplitude_summary(tls_amplitudes=tls_amps)
                tls_amp_dict.update({
                    'type' : 'alert',
                    'colour' : 'info',
                    'width' : 12,
                    'title' : 'All Amplitudes',
                    })

                # Get TLS amplitude summaries
                tls_ori_dict = self.create_tls_origins_summary(tls_origins=tls_coms)
                tls_ori_dict.update({
                    'type' : 'alert',
                    'colour' : 'info',
                    'width' : 12,
                    'title' : 'All TLS Origins',
                    })

                # Get images and format values
                image_dicts = []

                img = mf.get('level_uijs_pymol_by_group_png', {}).get(i_l+1, {}).get(n_g)
                if img:
                    image_dicts.append({
                        'width' : 4,
                        'title' : 'Average size (all datasets)',
                        'image' : self.image(img),
                        })

                img = mf.get('level_tls_amplitudes_png', {}).get(i_l+1, {}).get(n_g)
                if img:
                    image_dicts.append({
                        'width' : 4,
                        'title' : 'Amplitude distribution',
                        'image' : self.image(img),
                        })

                img = mf.get('level_uijs_pymol_by_group_rescaled_png', {}).get(i_l+1, {}).get(n_g)
                if img:
                    image_dicts.append({
                        'width' : 4,
                        'title' : 'Normalised (arbitrary scale)',
                        'image' : self.image(img),
                        })

                # Create panel dictionary
                group_tab = {
                    'id'        : level_tab['id']+'group{}'.format(self.counter.next()),
                    'alt_title' : 'Group {}'.format(n_g),
                    'title'     : 'Group {} - {}'.format(n_g, tls_obj.label),
                    'contents'  : [
                        {
                            'width':12,
                            'text':'<br>'.join(['Number of atoms: {}'.format(tls_obj.n_atoms)]),
                            },
                        ] + tls_mat_dicts_core + image_dicts + [tls_mat_dicts_dsts], #+ [tls_amp_dict, tls_ori_dict]
                    }
                tab_set['contents'].append(group_tab)

            tab_set['contents'][0].setdefault('active', True)

        return tabs

    def create_adp_level_tab(self):
        # -------------------------------->
        # Create tab for residual level
        # -------------------------------->
        residual_tab = {
            'id'        : 'lvlres',
            'alt_title' : 'Residual',
            'title'     : 'Final Level  (residual)',
            'contents'  : [],
            }
        tabs.append(residual_tab)
        # Get selection for fitted atoms
        atom_sel = flex.bool(p.atom_mask.tolist())
        # Create row for each residue
        for c_id in chain_ids:
            h = p.blank_master_hierarchy().select(atom_sel,copy_atoms=True)
            h = h.select(h.atom_selection_cache().selection('chain {}'.format(c_id)))

            # Panel for chain overview
            chain_image = fm.get_file('pml-residual-chain-template').format(c_id)
            stack_image = fm.get_file('png-residual-profile-template').format(c_id)
            aniso_image = fm.get_file('png-residual-anisotropy-template').format(c_id)
            panel = {
                'type'  : 'panel',
                'title' : 'Residual overview for chain {}'.format(c_id),
                'contents' : [
                    {
                        'width' : 8,
                        'title' : 'Fitted ADPs',
                        'image' : png2base64src_maybe(chain_image, print_on_missing=DEBUG),
                        'footnote' :'<strong>Please note</strong>: highly anisotropic atoms may not be visible in the above image.',
                        },
                    {
                        'width' : 6,
                        'title' : 'Fitted ADPs profile',
                        'image' : png2base64src_maybe(stack_image, print_on_missing=DEBUG),
                        },
                    {
                        'width' : 6,
                        'title' : 'Anisotropy by atom',
                        'image' : png2base64src_maybe(aniso_image, print_on_missing=DEBUG),
                        },
                    ],
                }
            residual_tab['contents'].append(panel)
            # Panel for each residue
            panel = {
                'type'  : 'panel',
                'title' : 'Residual components for chain {}'.format(c_id),
                'show'  : False,
                'contents' : [],
                }
            residual_tab['contents'].append(panel)
            for i_rg, rg in enumerate(h.residue_groups()):
                short_label = ShortLabeller.format(rg)
                long_label  = Labeller.format(rg)
                panel['contents'].append({'width':4, 'text':long_label})
                adp_image = fm.get_file('pml-residual-group-template').format(short_label)
                if os.path.exists(adp_image):
                    panel['contents'][-1]['image'] = self.image(adp_image)

    @classmethod
    def create_tls_matrix_summary(cls, tls_matrices):
        if tls_matrices.any():
            tls_mdl_str = ('<samp>' + \
                           'T: {T11:>9.3f}, {T22:>9.3f}, {T33:>9.3f}, {T12:>9.3f}, {T13:>9.3f}, {T23:>9.3f},<br>' + \
                           'L: {L11:>9.3f}, {L22:>9.3f}, {L33:>9.3f}, {L12:>9.3f}, {L13:>9.3f}, {L23:>9.3f},<br>' + \
                           'S: {S11:>9.3f}, {S12:>9.3f}, {S13:>9.3f}, {S21:>9.3f}, {S22:>9.3f}, {S23:>9.3f}, {S31:>9.3f}, {S32:>9.3f}, {S33:>9.3f}' + \
                           '</samp>').format(**tls_matrices.round(3)).replace(' ','&nbsp;')
        else:
            tls_mdl_str = "Zero-value TLS matrices"
        # Format output dictionary
        out_dict = {'type':'alert', 'colour':'info', 'text':tls_mdl_str, 'width':12}
        return out_dict

    @classmethod
    def create_tls_decomposition_summary(cls, tls_matrices, tolerance=1e-6, amplitude=1.0):
        T = tuple(tls_matrices[['T11','T22','T33','T12','T13','T23']])
        L = tuple(tls_matrices[['L11','L22','L33','L12','L13','L23']])
        S = tuple(tls_matrices[['S11','S12','S13','S21','S22','S23','S31','S32','S33']])
        a = amplitude
        assert a >= 0.0
        from mmtbx.tls.decompose import decompose_tls_matrices
        dcm = decompose_tls_matrices(T=T, L=L, S=S,
                                     l_and_s_in_degrees=True,
                                     tol=tolerance)
        #from IPython import embed; embed(); raise Exception()
        if dcm.is_valid():
            fmt_str_1 = '{:>8.6f}, {:>8.6f}, {:>8.6f}'
            dcm_lines = [
                    (
                        '<h5>Amplitudes of TLS Motions:</h5>'+\
                        '<samp>'+\
                        ('Vibration (&#8491;):  '    +fmt_str_1+'<br>').format(*[v*a for v in dcm.v_amplitudes])+\
                        ('Libration (&#176;):  '     +fmt_str_1+'<br>').format(*[cls.RAD2DEG*v*a for v in dcm.l_amplitudes])+\
                        ('Screw   (&#8491;/&#176;):  ' +fmt_str_1+'<br>').format(*[v/cls.RAD2DEG for v in dcm.s_amplitudes])+\
                        '<br>'+\
                        ('Libration (Rad):  '          +fmt_str_1+'<br>').format(*[v*a for v in dcm.l_amplitudes])+\
                        ('Screw   (&#8491;/Rad):  '    +fmt_str_1+'<br>').format(*[v*a for v in dcm.s_amplitudes])+\
                        '</samp>'
                    ),
                    (
                        '<h5>Directions/Axes of TLS Motions:</h5>'+\
                        '<samp>'+\
                        (
                            'Vibration Axes:  ({:>9.3f}, {:>9.3f}, {:>9.3f})<br>'+\
                            '                 ({:>9.3f}, {:>9.3f}, {:>9.3f})<br>'+\
                            '                 ({:>9.3f}, {:>9.3f}, {:>9.3f})<br>'
                        ).format(*numpy.concatenate(dcm.v_axis_directions))+'<br>'+\
                        (
                            'Libration Axes:  ({:>9.3f}, {:>9.3f}, {:>9.3f})<br>'+\
                            '                 ({:>9.3f}, {:>9.3f}, {:>9.3f})<br>'+\
                            '                 ({:>9.3f}, {:>9.3f}, {:>9.3f})<br>'
                        ).format(*numpy.concatenate(dcm.l_axis_directions))+'<br>'+\
                        (
                            'Libration Axis   ({:>9.3f}, {:>9.3f}, {:>9.3f})<br>'+\
                            'Intersections:   ({:>9.3f}, {:>9.3f}, {:>9.3f})<br>'+\
                            '                 ({:>9.3f}, {:>9.3f}, {:>9.3f})<br>'
                        ).format(*numpy.concatenate(dcm.l_axis_intersections)) +\
                        '</samp>'
                    ),
                        ]
            dcm_string = ('<hr>'.join(dcm_lines)).replace(' ','&nbsp;')
            colour = 'success'
        else:
            dcm_string = 'Unable to decompose TLS matrices: <br>&nbsp;&nbsp;&nbsp;&nbsp;{}'.format(dcm.error())
            colour = 'danger'
        # Format output dictionary
        out_dict = {'type':'alert', 'colour':colour, 'text':dcm_string, 'width':12}
        return out_dict

    @classmethod
    def create_tls_amplitude_summary(cls, tls_amplitudes):
        tls_amp_table = tls_amplitudes.T
        # Change the names of the "Mode" columns
        new_columns = ['Mode {}'.format(m) for m in tls_amp_table.columns]
        tls_amp_table.columns = new_columns
        # Remove columns names and add index name
        tls_amp_table.columns.name = 'Dataset'
        html_table = tls_amp_table.to_html(
                float_format=lambda v: '{:.6f}'.format(v),
                bold_rows=False,
                classes=['table table-condensed table-hover datatable nowrap text-center'],
                border=0)
        html_table = html_table.replace('<th>', '<th class="text-center">')
        return {'table': html_table}

    @classmethod
    def create_tls_origins_summary(cls, tls_origins):
        # Change the names of the columns
        tls_origins.columns = ['x','y','z']
        # Remove columns names and add index name
        tls_origins.columns.name = 'Dataset'
        html_table = tls_origins.to_html(
                float_format=lambda v: '{:.3f}'.format(v),
                bold_rows=False,
                classes=['table table-condensed table-hover datatable nowrap text-center'],
                border=0)
        html_table = html_table.replace('<th>', '<th class="text-center">')
        return {'table': html_table}
