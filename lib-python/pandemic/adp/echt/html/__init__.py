from libtbx import adopt_init_args

import numpy, pandas

from pandemic.adp.html import (
    divs,
    HtmlSummary,
    )

from pandemic.adp import constants


class EchtModelHtmlSummary(HtmlSummary):


    def __init__(self,
        hierarchy_summary_task,
        model_summary_task,
        model_object,
        isotropic_mask,
        parameters,
        ):
        adopt_init_args(self, locals())

    def short_summary(self):

        hf = self.hierarchy_summary_task.result.output_files
        mf = self.model_summary_task.result.output_files

        img1 = mf.get('all_levels_uijs_profiles_png', {})
        img2 = mf.get('all_levels_uijs_anisotropy_png', {})
        img3 = hf.get('level_partitions_png', {})

        chain_ids = sorted(set(list(img1.keys())+list(img2.keys())+list(img3.keys())))

        b_factor_table_blocks = self.format_b_factor_table_blocks(chains=None)

        output = divs.Panel(title='Output: Chain-by-Chain Disorder Summaries')

        ###

        block = divs.Block(
            contents = self.format_summary(
                (
                    'Visualise in pymol by running: ' +
                    self.wrap_string('pymol {}', 'pre')
                ).format(mf.get("pymol_script")),
                type = 'none',
            ),
        )
        output.append(block)

        ###

        for c_id in chain_ids:
            block = divs.Block(
                title = 'Chain {}'.format(c_id),
                contents = [
                    divs.Block(
                        width = 4,
                        text = 'Level-by-level ADP Profile',
                        image = self.image(img1.get(c_id)),
                    ),
                    divs.Block(
                        width = 4,
                        text = 'Total model anisotropy',
                        image = self.image(img2.get(c_id)),
                    ),
                    divs.Block(
                        width = 4,
                        text = 'Partition Schematic',
                        image = self.image(img3.get(c_id)),
                    ),
                ],
            )
            output.append(block)

        return b_factor_table_blocks + [output]

    def main_summary(self):

        mf = self.model_summary_task.result.output_files

        output = divs.Tab(
            id = 'levels',
            title = 'ECHT level-by-level TLS parameterisation',
            alt_title = 'Hierarchical Model Summary',
        )

        output.extend(
            self.model_hierarchy_summary_string_list(
                model_object = self.model_object,
                isotropic_mask = self.isotropic_mask,
            )
        )
        output.extend(
            self.model_parameter_summary_string_list(
                model_object = self.model_object,
                isotropic_mask = self.isotropic_mask,
            )
        )

        output.append(
            divs.Alert(
                title = 'Level B-factor Distributions (all chains)',
                width = 12,
                image = self.image(
                    mf.get('b_factor_distributions', None)
                    ),
                )
            )

        output.extend(
            self.format_b_factor_table_blocks(chains='all', split_tables=True)
        )

        overview_tab = self.overview_tab(parent_tab_id=output.id)
        tls_level_tabs = self.create_tls_level_tabs(parent_tab_id=output.id)
        adp_level_tab = self.create_adp_level_tab(parent_tab_id=output.id)

        tab_set = divs.TabSet(
            contents = [overview_tab] + tls_level_tabs + [adp_level_tab],
        )
        tab_set.set_active()

        output.append(tab_set)

        return [output]

    def model_hierarchy_summary_string_list(self,
        model_object,
        isotropic_mask,
        ):

        mo = model_object

        n_tls_groups = numpy.sum([len(gs) for gs in mo.tls_objects])

        atoms_per_group_per_level = [[g.n_atoms for g in gs] for gs in mo.tls_objects]
        group_sizes_per_level = ['{} - {} atoms'.format(min(n), max(n)) for n in atoms_per_group_per_level]

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
            s += '\nAnisotropic - Isotropic atoms: {} - {}'.format(isotropic_mask.n_anisotropic, isotropic_mask.n_isotropic)
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
        else:
            n_iso = isotropic_mask.n_isotropic
            n_ani = isotropic_mask.n_anisotropic

        n_adp_params = 6 * n_ani + n_iso
        n_tls_params = sum(n_tls_params_per_level)
        n_params = n_adp_params + n_tls_params

        s = ''
        s += '\n> Parameterisation Summary:'
        s += '\nParameters per group: {} (20 from matrices + {} from amplitudes)'.format(n_tls_params_per_group, mo.n_datasets-1)
        s += '\nParameters per TLS level:'
        for i_l, l in enumerate(mo.tls_level_names):
            s += '\n\tLevel {} ({}): {}'.format(i_l+1, l, n_tls_params_per_level[i_l])
        s += '\nTotal TLS parameters: {}'.format(n_tls_params)
        if adp_level_opt:
            s += '\nTotal ADP ({} level) parameters: {} (6*{} + 1*{})'.format(mo.adp_level_name, n_adp_params, n_ani, n_iso)
        s += '\nTotal parameters: {} ({:.2f} per atom per dataset)'.format(n_params, float(n_params)/float(mo.n_atoms * mo.n_datasets))

        return self.format_summary(s, width=6)

    def overview_tab(self, parent_tab_id):
        """Make overview tab for the hierarchical model"""

        hf = self.hierarchy_summary_task.result.output_files
        mf = self.model_summary_task.result.output_files

        img1 = mf.get('all_levels_uijs_profiles_png', {})
        img2 = mf.get('all_levels_uijs_anisotropy_png', {})

        img3 = mf.get('level_uijs_pymol_by_chain_png', {})
        img4 = mf.get('level_uijs_profiles_png', {})
        img5 = mf.get('level_uijs_anisotropy_png', {})

        chain_ids = sorted(img1.keys())

        # -------------------------------->
        # Create overview sub-tab
        # -------------------------------->
        overview_tab = divs.Tab(
            id = parent_tab_id+'overview',
            title = 'Overview of the parameterised hierarchical ADP model',
            alt_title = 'Overview',
        )

        ###

        block = divs.Block(
            contents = self.format_summary(
                (
                    'Visualise in pymol by running: ' +
                    self.wrap_string('pymol {}', 'pre')
                ).format(mf.get("pymol_script")),
                type='none',
            ),
        )
        overview_tab.append(block)

        ###

        tab_set = overview_tab.append(divs.TabSet())

        # Split the panels up by chain
        for c_id in chain_ids:

            # Split up the chains with divider panels
            tab = divs.Tab(
                title = 'Levels for Chain {}'.format(c_id),
                alt_title = 'Chain {}'.format(c_id),
            )
            tab_set.append(tab)

            blocks = self.format_b_factor_table_blocks(chains=[c_id], split_tables=True)
            tab.extend(blocks)

            block = divs.Alert(
                contents = [
                    divs.Block(
                        width = 6,
                        title = 'Level-by-level Disorder Profile',
                        image = self.image(img1.get(c_id)),
                    ),
                    divs.Block(
                        width = 6,
                        title = 'Anisotropy of Total Disorder',
                        image = self.image(img2.get(c_id)),
                    ),
                ],
            )
            tab.append(block)

            n_levels = self.model_object.n_levels
            for i_l, l in enumerate(self.model_object.all_level_names):

                block = divs.Alert(
                    title = 'Level {} of {} ({})'.format(i_l+1, n_levels, l),
                    contents = [
                        divs.Block(
                            width = 4,
                            image = self.image(img3.get(l, {}).get(c_id)),
                        ),
                        divs.Block(
                            width = 4,
                            image = self.image(img4.get(l, {}).get(c_id)),
                        ),
                        divs.Block(
                            width = 4,
                            image = self.image(img5.get(l, {}).get(c_id)),
                        ),
                    ],
                )
                tab.append(block)

        tab_set.set_active()

        return overview_tab

    def create_tls_level_tabs(self, parent_tab_id):

        hf = self.hierarchy_summary_task.result.output_files
        mf = self.model_summary_task.result.output_files

        mo = self.model_object

        img1 = hf.get('level_atoms_pymol_png', {})
        img2 = mf.get('level_uijs_pymol_by_chain_png', {})
        img3 = mf.get('level_uijs_profiles_png', {})
        img4 = mf.get('level_uijs_anisotropy_png', {})

        chain_ids = []
        if list(img1.values()):
            chain_ids += list(img1.values())[0].keys()
        if list(img2.values()):
            chain_ids += list(img2.values())[0].keys()
        if list(img3.values()):
            chain_ids += list(img3.values())[0].keys()
        if list(img4.values()):
            chain_ids += list(img4.values())[0].keys()
        chain_ids = sorted(set(chain_ids))

        # -------------------------------->
        # Create tab for each level
        # -------------------------------->
        tab_list = []
        for i_l, l in enumerate(mo.tls_level_names):

            # Create dictionary for this tab and add to tab_list
            level_tab = divs.Tab(
                id = parent_tab_id+'lvl{}'.format(i_l),
                alt_title = 'Level {} ({})'.format(i_l+1, l),
            )
            tab_list.append(level_tab)

            txt = """
            > Level {l_num} ({l_name})
            Level {l_num} of {n_levels}. Composed of {n_groups} groups.
            """.format(
                l_num = i_l+1,
                l_name = l,
                n_levels = mo.n_levels,
                n_groups = len(mo.tls_objects[i_l]),
            )
            
            level_tab.extend(
                self.format_summary(txt, classes=['square-corners-top'])
                )

            level_tab.append(
                divs.Alert(
                    title = "Level B-factor distribution",
                    width = 12,
                    image = self.image(
                        mf.get('b_factor_distributions_level',{}).get(l)
                        ),
                    )
                )

            ########################################################

            # Add overview at the top of the tab
            tab_set = divs.TabSet(title="Chain-by-Chain Summaries")
            level_tab.append(tab_set)

            # Make tabs
            for c_id in chain_ids:

                chain_tab = divs.Tab(
                    id = level_tab['id']+'chain{}'.format(self.counter.next()),
                    title = 'Disorder Summary for Chain {}'.format(c_id),
                    alt_title = 'Chain {}'.format(c_id),
                    contents = [
                        divs.Block(
                            width = 6,
                            title = 'Image of chain, coloured by group',
                            image = self.image(img1.get(l,{}).get(c_id)),
                        ),
                        divs.Block(
                            width = 6,
                            title = 'Level ADPs' + ' (average size over all datasets)'*(mo.n_datasets>1),
                            image = self.image(img2.get(l,{}).get(c_id)),
                            footnote = '<strong>Please note</strong>: highly anisotropic atoms may not be visible in the above image.',
                        ),
                        divs.Block(
                            width = 6,
                            title = 'Level ADPs by mode' + ' (average size over all datasets)'*(mo.n_datasets>1),
                            image = self.image(img3.get(l,{}).get(c_id)),
                        ),
                        divs.Block(
                            width = 6,
                            title = 'ADP Anisotropy' + ' (from average over all datasets)'*(mo.n_datasets>1),
                            image = self.image(img4.get(l,{}).get(c_id)),
                        ),
                    ],
                )
                tab_set.append(chain_tab)
            tab_set.set_active()

            ########################################################

            # Read in the TLS modes and amplitudes for this level
            tls_matrices   = pandas.read_csv(mf['level_tls_matrices_csv'][l]).set_index(['group','mode']).drop('Unnamed: 0', axis=1, errors='ignore')
            tls_amplitudes = pandas.read_csv(mf['level_tls_amplitudes_csv'][l]).set_index(['group','mode']).drop('Unnamed: 0', axis=1, errors='ignore')

            ########################################################

            # Amplitude summaries
            tls_amp_panel = divs.Panel(
                title = 'Group Amplitudes',
                contents = [
                    # Note about the scaling of the amplitudes in the table
                    divs.Block(text='Amplitudes are in the scale of B-factors'),
                ],
            )
            level_tab.append(tls_amp_panel)

            # Get TLS amplitude summaries
            tls_amplitudes_table_data = tls_amplitudes.copy()
            av_label = 'average'
            tls_amplitudes_table_data[av_label] = tls_amplitudes_table_data[mo.dataset_labels].mean(axis=1)
            # Reorder columns & multiply all of the flot columns by EIGHTPISQ
            cols = tls_amplitudes_table_data.columns.tolist()
            for c in [av_label] + mo.dataset_labels:
                # remove from list so can be reordered later
                cols.remove(c)
                # scale column in table
                tls_amplitudes_table_data[c] *= constants.EIGHTPISQ
            # resort columns
            tls_amplitudes_table_data = tls_amplitudes_table_data[cols+[av_label]+mo.dataset_labels]

            # Create html object
            tls_amp_data = self.create_tls_amplitude_summary(tls_amplitudes=tls_amplitudes_table_data)
            tls_amp_panel.append(tls_amp_data)

            ########################################################

            tab_set = divs.TabSet(
                title = "Group-by-Group Summaries",
                classes = ['nav-hover'],
            )
            level_tab.append(tab_set)

            # Extract groups for each level
            for i_g, tls_obj in enumerate(mo.tls_objects[i_l]):
                n_g = i_g + 1

                # Create panel dictionary
                group_tab = divs.Tab(
                    id = level_tab['id']+'group{}'.format(self.counter.next()),
                    alt_title = 'Group {}'.format(n_g),
                )
                tab_set.append(group_tab)

                txt = """
                > Group {g_num} - {g_name}
                Number of atoms: {n_atoms}
                """.format(
                    g_num = n_g,
                    g_name = tls_obj.label,
                    n_atoms = tls_obj.n_atoms,
                )
                group_tab.extend(self.format_summary(txt, classes=['square-corners-top']))

                # Get images and format values
                image_block = divs.Block()

                img = mf.get('level_uijs_pymol_by_group_png', {}).get(i_l+1, {}).get(n_g)
                if img:
                    image_block.append(
                        divs.Block(
                            width = 4,
                            title = 'Average size (all datasets)',
                            image = self.image(img),
                        )
                    )

                img = mf.get('level_tls_amplitudes_png', {}).get(i_l+1, {}).get(n_g)
                if img:
                    image_block.append(
                        divs.Block(
                            width = 4,
                            title = 'Amplitude distribution',
                            image = self.image(img),
                        )
                    )

                img = mf.get('level_uijs_pymol_by_group_rescaled_png', {}).get(i_l+1, {}).get(n_g)
                if img:
                    image_block.append(
                        divs.Block(
                            width = 4,
                            title = 'Normalised (arbitrary scale)',
                            image = self.image(img),
                        )
                    )

                if image_block.contents:
                    group_tab.append(image_block)

                # Extract TLS values for this group
                tls_mats = tls_matrices.loc[n_g].drop(columns=['label'])
                tls_amps = tls_amplitudes.loc[n_g].drop(columns=['label'])
                tls_amps_av = tls_amps.mean(axis=1).values # average over datasets
                tls_coms = pandas.DataFrame(index=tls_amps.T.index, data=numpy.array(tls_obj.origins))

                # Div for the TLS matrices summary
                tls_mat_block = divs.ScrollX(classes=['bordered'])
                group_tab.append(tls_mat_block)

                # Div for the average over all the datasets
                core_block = divs.Block()
                tls_mat_block.append(core_block)
                # Header
                core_block.append(
                    divs.Block(
                        title = 'TLS matrices averaged over all datasets',
                        classes = ['text-center'],
                        text = (
                            ' <strong>scroll for individual datasets</strong> '
                            '<i class="fa fa-angle-double-right fa-2x" style="vertical-align: middle;"></i>'
                            '<br>'
                        ),
                    )
                )
                # Average amplitudes for each of the modes
                core_block.append(
                    divs.Alert(
                        width = 6,
                        title = 'TLS amplitudes',
                        text = '<br>'.join([
                            '<samp>Mode {}: <strong>{:>.5f} &#8491;&#178; ({:>.3f} B-factor)</strong></samp>'.format(
                                i+1,
                                v,
                                constants.EIGHTPISQ*v,
                            ) for i,v in enumerate(tls_amps_av)
                        ]),
                    )
                )
                # TLS Origin - NONE HERE
                core_block.append(
                    divs.Alert(
                        width = 6,
                        title = 'TLS Group Origin',
                        text = 'No origin has been applied to the libration axes',
                    )
                )
                # Matrices & Decompositions
                for i_mode, (idx, mats) in enumerate(tls_mats.iterrows()):
                    # Average amplitude for the mode
                    mode_amp_av = tls_amps_av[i_mode]
                    # Matrix summary
                    block = self.create_tls_matrix_summary(tls_matrices=mode_amp_av*mats)
                    block.width = 12
                    block.title = 'Mode {} (average over all datasets)'.format(idx)
                    core_block.append(block)
                    # TLS decomposition summary
                    block = self.create_tls_decomposition_summary(
                        tls_matrices = mats,
                        tls_origin = None,
                        tls_amplitude = mode_amp_av,
                        tolerance = self.parameters.model.echt.tolerances.tls_matrix_tolerance,
                        )
                    block.width = 12
                    block.title = 'Average TLS decomposition of mode {}:'.format(i_mode+1)
                    core_block.append(block)

                # Now add the dicts for the datasets
                for i_dst, (dst_lab, dst_amps) in enumerate(tls_amps.T.iterrows()):

                    orgn = tuple(tls_coms.loc[dst_lab].values)

                    # Div for the each dataset
                    dset_block = divs.Block()
                    tls_mat_block.append(dset_block)
                    # Header
                    dset_block.append(
                        divs.Block(
                            title = 'TLS matrices for dataset "{}"'.format(dst_lab),
                            classes = ['text-center'],
                            text = (
                                '<i class="fa fa-angle-double-left fa-2x" style="vertical-align: middle;"></i>'
                                ' <strong>scroll for other datasets</strong> '
                                '<i class="fa fa-angle-double-right fa-2x" style="vertical-align: middle;"></i>'
                                '<br>'
                            ),
                        )
                    )
                    # Amplitudes for each of the modes
                    dset_block.append(
                        divs.Alert(
                            width = 6,
                            title = 'TLS amplitudes',
                            text = '<br>'.join([
                                '<samp>Mode {}: <strong>{:>.5f} &#8491;&#178; ({:>.3f} B-factor)</strong></samp>'.format(
                                    i+1,
                                    v,
                                    constants.EIGHTPISQ*v,
                                ) for i,v in enumerate(dst_amps.values)
                            ]),
                        )
                    )
                    # TLS Origin
                    dset_block.append(
                        divs.Alert(
                            width = 6,
                            title = 'TLS Group Origin',
                            text = '<samp>({:>.3f}, {:>.3f}, {:>.3f})</samp>'.format(*orgn),
                        )
                    )
                    # Matrices & Decompositions
                    for i_mode, (idx, mats) in enumerate(tls_mats.iterrows()):
                        # TLS matrix summary
                        block = self.create_tls_matrix_summary(tls_matrices=dst_amps.values[i_mode]*mats)
                        block.width = 12
                        block.title = 'Mode {} (with amplitudes):'.format(i_mode+1)
                        dset_block.append(block)
                        # TLS decomposition summary
                        block = self.create_tls_decomposition_summary(
                            tls_matrices = mats,
                            tls_origin = orgn,
                            tls_amplitude = dst_amps.values[i_mode],
                            tolerance = self.parameters.model.echt.tolerances.tls_matrix_tolerance,
                            )
                        block.width = 12
                        block.title = 'TLS decomposition of mode {}:'.format(i_mode+1)
                        dset_block.append(block)

                # Update the last block
                last_block_top = tls_mat_block.contents[-1].contents[0]
                last_block_top.text = last_block_top.text.replace('<i class="fa fa-angle-double-right fa-2x" style="vertical-align: middle;"></i>','')

            tab_set.set_active()

        return tab_list

    def create_adp_level_tab(self, parent_tab_id):
        """Create tab for residual level"""

        hf = self.hierarchy_summary_task.result.output_files
        mf = self.model_summary_task.result.output_files

        mo = self.model_object

        l = mo.adp_level_name
        i_l = mo.all_level_names.index(l)

        img1 = mf.get('level_uijs_pymol_by_chain_png', {}).get(l, {})
        img2 = mf.get('level_uijs_profiles_png', {}).get(l, {})
        img3 = mf.get('level_uijs_anisotropy_png', {}).get(l, {})

        chain_ids = list(img1.keys())+list(img2.keys())+list(img3.keys())
        chain_ids = sorted(set(chain_ids))

        atomic_tab = divs.Tab(
            id = parent_tab_id+'lvlatm',
            alt_title = 'Level {} ({})'.format(i_l+1, l),
        )

        txt = """
        > Level {l_num} ({l_name})
        Level {l_num} of {n_levels}. Composed of {n_atoms} atoms.
        """.format(
            l_num = i_l+1,
            l_name = l,
            n_levels = mo.n_levels,
            n_atoms = mo.n_atoms,
        )

        atomic_tab.extend(
            self.format_summary(txt, classes=['square-corners-top'])
            )

        atomic_tab.append(
            divs.Alert(
                title = "Level B-factor distribution",
                width = 12,
                image = self.image(
                    mf.get('b_factor_distributions_level',{}).get(l)
                    ),
                )
            )

        # Add overview at the top of the tab
        tab_set = divs.TabSet(title="Chain-by-Chain Summaries")
        atomic_tab.append(tab_set)

        for c_id in chain_ids:

            chain_tab = divs.Tab(
                id = atomic_tab['id']+'chain{}'.format(self.counter.next()),
                title = 'Summary for Chain {}'.format(c_id),
                alt_title = 'Chain {}'.format(c_id),
            )
            tab_set.append(chain_tab)

            chain_tab.append(
                divs.Block(
                    contents = [
                        divs.Block(width=3),
                        divs.Block(
                            width = 6,
                            title = 'Level ADPs',
                            image = self.image(img1.get(c_id)),
                            footnote = '<strong>Please note</strong>: highly anisotropic atoms may not be visible in the above image.',
                        ),
                    ],
                )
            )

            chain_tab.append(
                divs.Block(
                    contents = [
                        divs.Block(
                            width = 6,
                            title = 'Level ADPs Profile',
                            image = self.image(img2.get(c_id)),
                        ),
                        divs.Block(
                            width = 6,
                            title = 'ADP Anisotropy Profile',
                            image = self.image(img3.get(c_id)),
                        ),
                    ],
                )
            )

        tab_set.set_active()

        return atomic_tab

    @classmethod
    def create_tls_matrix_summary(cls,
        tls_matrices,
        ):

        if tls_matrices.any():

            tls_mdl_str = (
                '<strong>T Matrix:</strong><br><br>' + \
                '{T11:>7.3f}, {T12:>7.3f}, {T13:>7.3f},<br>' + \
                '{BLANKSP:^8} {T22:>7.3f}, {T23:>7.3f},<br>' + \
                '{BLANKSP:^8} {BLANKSP:^8} {T33:>7.3f},<br><br><br>' + \
                '<strong>L Matrix:</strong><br><br>' + \
                '{L11:>7.3f}, {L12:>7.3f}, {L13:>7.3f},<br>' + \
                '{BLANKSP:^8} {L22:>7.3f}, {L23:>7.3f},<br>' + \
                '{BLANKSP:^8} {BLANKSP:^8} {L33:>7.3f},<br><br><br>' + \
                '<strong>S Matrix:</strong><br><br>' + \
                '{S11:>7.3f}, {S12:>7.3f}, {S13:>7.3f},<br>' + \
                '{S21:>7.3f}, {S22:>7.3f}, {S23:>7.3f},<br>' + \
                '{S31:>7.3f}, {S32:>7.3f}, {S33:>7.3f}'
            ).format(
                BLANKSP='-',
                **tls_matrices.round(3)
            ).replace(' ','&nbsp;').split('<br><br><br>')
            block_width = 4

        else:

            tls_mdl_str = ["Zero-value TLS matrices"]
            block_width = 12

        # Format output dictionary
        out_block = divs.Alert(width=12)

        for b in tls_mdl_str:
            out_block.append(
                divs.Block(
                    width = block_width,
                    text = '<samp><strong>'+b+'</strong></samp>',
                    classes = ['bordered'],
                )
            )

        return out_block

    @classmethod
    def create_tls_decomposition_summary(cls,
        tls_matrices,
        tls_origin = None,
        tls_amplitude = 1.0,
        tolerance = None,
        ):

        T = tuple(tls_matrices[['T11','T22','T33','T12','T13','T23']])
        L = tuple(tls_matrices[['L11','L22','L33','L12','L13','L23']])
        S = tuple(tls_matrices[['S11','S12','S13','S21','S22','S23','S31','S32','S33']])

        # Actual multiplier is the square root of the amplitude. Used
        # on the vibrational and librational amplitudes ONLY.
        # Screw component amplitudes is unaffected by multipliers!
        assert tls_amplitude >= 0.0
        multiplier = (tls_amplitude ** 0.5)

        if tls_origin is None:
            tls_origin = numpy.zeros(3)
        tls_origin = numpy.array(tls_origin)
        assert len(tls_origin) == 3

        from mmtbx.tls.decompose import decompose_tls_matrices
        dcm = decompose_tls_matrices(
            T = T,
            L = L,
            S = S,
            l_and_s_in_degrees = True,
            tol = tolerance,
            )

        if dcm.is_valid():

            fmt_str_1 = '{:>5.3f}'
            fmt_str_2 = '{0:>5.2f}, {1:>5.2f}, {2:>5.2f}'

            vibration_block = ""

            for i in range(2,-1,-1):
                va_format = fmt_str_1.format(multiplier*dcm.v_amplitudes[i])
                vd_format = fmt_str_2.format(*dcm.v_axis_directions[i])
                vibration_block += 'Vibration Amplitude: <strong>{a} &#8491;</strong><br>'.format(a=va_format)
                vibration_block += 'Axis Direction: <strong>({d})</strong><br><br>'.format(d=vd_format)

            libration_block = ""

            for i in range(2,-1,-1):
                la_format = fmt_str_1.format(multiplier*constants.RAD2DEG*dcm.l_amplitudes[i])
                sa_format = fmt_str_1.format(dcm.s_amplitudes[i]/constants.RAD2DEG)
                ld_format = fmt_str_2.format(*dcm.l_axis_directions[i])
                li_format = fmt_str_2.format(*(tls_origin+dcm.l_axis_intersections[i]))

                libration_block += 'Libration Amplitude: <strong>{a}&#176;</strong><br>'.format(a=la_format)
                libration_block += 'Screw Amplitude: <strong>{a} &#8491;/&#176;</strong><br>'.format(a=sa_format)
                libration_block += 'Axis Direction:    <strong>({d})</strong><br>'.format(d=ld_format)
                libration_block += 'Axis Intersection: <strong>({d})</strong><br><br>'.format(d=li_format)

            vibration_block = ("<samp>\n" + vibration_block + "\n</samp>").replace(' ','&nbsp;')
            libration_block = ("<samp>\n" + libration_block + "\n</samp>").replace(' ','&nbsp;')

            dcm_blocks = [
                divs.Block(
                    width = 6,
                    title = "Vibrational Modes",
                    text = vibration_block,
                ),
                divs.Block(
                    width = 6,
                    title = "Librational &amp; Screw Modes",
                    text = libration_block,
                ),
            ]
            colour = None

        else:

            dcm_blocks = [
                divs.Block(
                    text = 'Unable to decompose TLS matrices: <br>&nbsp;&nbsp;&nbsp;&nbsp;{}'.format(dcm.error()).replace(' ','&nbsp;'),
                ),
            ]
            colour = 'danger'

        # Format output dictionary
        out_block = divs.Alert(
            title = 'Descriptions of TLS Motions',
            contents = dcm_blocks,
        )
        if colour is not None:
            out_block.colour = colour

        return out_block

    @classmethod
    def create_tls_amplitude_summary(cls, tls_amplitudes):
        html_table = tls_amplitudes.reset_index().to_html(
            index=False,
            float_format=lambda v: '{:.3f}'.format(v),
            bold_rows=False,
            classes=['table table-hover datatable nowrap text-center'],
            justify='center',
            border=0,
            )
        html_table = html_table \
            .replace('<th>', '<th class="text-center">') \
            .replace('border="1" ', '')
        return divs.Block(table=html_table)

    @classmethod
    def create_tls_origins_summary(cls, tls_origins):
        # Change the names of the columns
        tls_origins.columns = ['x','y','z']
        # Remove columns names and add index name
        tls_origins.columns.name = 'Dataset'
        html_table = tls_origins.to_html(
            float_format=lambda v: '{:.3f}'.format(v),
            bold_rows=False,
            classes=['table table-hover datatable nowrap text-center'],
            justify='center',
            border=0,
            )
        html_table = html_table \
            .replace('<th>', '<th class="text-center">') \
            .replace('border="1" ', '')
        return divs.Block(table=html_table)

    def format_b_factor_table_blocks(self, chains=None, split_tables=False):

        b_table = self.model_summary_task.result.level_b_factor_statistics_table

        levels_ordered = b_table[b_table['Chain'] == 'all']['Level'].tolist()

        # Which columns to use in the output table
        output_cols = [c for c in b_table.columns if c not in ['Level', 'Chain']]

        if (chains is not None):
            if isinstance(chains, str):
                chains = [chains]
            table_sel = b_table['Chain'].isin(chains)
            b_table = b_table[table_sel]

        output_blocks = []

        if (split_tables is False):
            output_val_sets = [output_cols]
        else:
            output_val_sets = [[c] for c in output_cols]

        for val_set in output_val_sets:

            # Pivot table
            table = b_table.pivot(index='Level', columns='Chain', values=val_set)

            # Don't need chain level if only one chain
            if (chains is not None) and (len(chains) == 1):
                table = table.droplevel('Chain', axis=1)

            # Return levels to correct order (sorted alphabetically after pivoting)
            table = table.reindex(levels_ordered)

            title = 'Average B-factors for Levels'
            if chains == ['all']:
                title += ' (all chains)'
            elif (chains is not None):
                title += ' (chain{s} {ids})'.format(
                    s = 's' if len(chains) > 1 else '',
                    ids = ' & '.join(chains),
                )

            table_block = divs.Alert(
                title = title,
                width = 12 if (split_tables is False) else 6,
                table = table.round(1)\
                    .to_html(index=True, bold_rows=False, classes=['table table-hover nowrap'])\
                    .replace('border="1" ', ''),
            )

            output_blocks.append(table_block)

        return output_blocks
