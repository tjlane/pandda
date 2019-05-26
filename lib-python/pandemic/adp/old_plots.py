    # adp.analysis
    def alpha_scatter(self,
        prefix,
        hierarchies,
        title=None,
        y_lab='Isotropic B ($\AA^2$)',
        y_lim=None,
        v_line_hierarchy=None,
        rotate_x_labels=True,
        plot_mean=True,
        add_jitter=True,
        ):
        """Plot scatter plot over atoms (by residue) with alpha to prevent overplotting for a series of hierarchies (plotted values are the average B-factors of each residue of the hierarchies)"""

        m_h = hierarchies[0]

        # Check all hierarchies are the same
        for h in hierarchies:
            assert m_h.is_similar_hierarchy(h)

        colours = pyplot.cm.rainbow(numpy.linspace(0,1,len(hierarchies)))

        # Create a plot for each chain
        for chain_id in sorted(set([c.id for c in m_h.chains()])):

            # Create a selection for each chain and select it in each hierarchy
            sel = m_h.atom_selection_cache().selection('chain {}'.format(chain_id))
            sel_hs = [h.select(sel) for h in hierarchies]
            sel_mh = sel_hs[0]

            # Filename!
            filename = prefix + '-chain_{}.png'.format(chain_id)

            # Create x-values for each residue starting from 1
            x_vals = numpy.array(range(len(list(sel_mh.residue_groups()))))+1
            x_labels = ['']+[ShortLabeller.format(rg) for rg in sel_mh.residue_groups()]

            # Create the output figure
            fig, axis = pyplot.subplots(nrows=1, ncols=1)
            if title is not None: axis.set_title(label=str(title))
            axis.set_xlabel('Residue')
            axis.set_ylabel(y_lab)
            if y_lim:
                axis.set_ylim(y_lim)

            # Alpha value for scatter -- TODO optimise value?
            alpha = 1.0 / (1.0 + float(len(hierarchies))**0.5)

            # Mean values to be plotted as line
            mean_arr = numpy.zeros(x_vals.shape)

            # Iterative though hierarchies
            for i_h, h in enumerate(sel_hs):
                # Extract b-factors from this hierarchy
                y_vals = [numpy.array(rg.atoms().extract_b()) for rg in h.residue_groups()]
                # Append to mean arr
                y_mean = numpy.array([y.mean() for y in y_vals])
                assert len(y_mean) == len(mean_arr)
                mean_arr = mean_arr + y_mean
                # Reformat x,y for plotting
                y = numpy.concatenate(y_vals)
                x = numpy.concatenate([[x]*len(yvs) for x, yvs in zip(x_vals, y_vals)])
                if add_jitter is True:
                    jitter = numpy.random.randn(len(x)) * 0.1 / 3.0 # Normalise to be approx width of bar
                    jitter = jitter - jitter.mean()
                else:
                    jitter = 0.0
                # Plot the scatter
                axis.scatter(x=x+jitter, y=y, s=5, alpha=alpha, facecolor=colours[i_h], lw=0.0)
                # Plot line
                axis.plot(x_vals, [v.mean() for v in y_vals], '-', color=colours[i_h], alpha=0.75, zorder=5, lw=0.75)

            # Skip if no values plotted
            if (mean_arr == 0.0).all():
                continue

            if plot_mean is True:
                mean_arr = mean_arr / float(len(hierarchies))
                axis.plot(x_vals, mean_arr, 'k-', lw=1.0, zorder=6)

            if v_line_hierarchy is not None:
                v_lines = numpy.where(numpy.array([max(rg.atoms().extract_b()) for rg in v_line_hierarchy.select(sel).residue_groups()], dtype=bool))[0] + 1.5
                for val in v_lines:
                    axis.axvline(x=val, c='grey', ls='solid', label='boundaries', lw=0.5, zorder=1)
                # Add a line at zero
                h = axis.axvline(x=0.5, c='grey', ls='solid', label='boundaries', lw=0.5, zorder=1)

            # Axis ticks & labels
            x_ticks = numpy.arange(1, len(x_labels)+1, int(max(1.0, numpy.floor(len(x_labels)/20))))
            #x_ticks = numpy.unique(numpy.floor(numpy.linspace(1,len(x_labels)+1,20)))
            axis.set_xticks(x_ticks)
            axis.set_xticklabels([x_labels[int(i)] if (i<len(x_labels)) and (float(int(i))==i) else '' for i in axis.get_xticks()])
            # Rotate axis labels
            if rotate_x_labels: pyplot.setp(axis.get_xticklabels(), rotation=90)

            # Format and save
            fig.tight_layout()
            fig.savefig(filename,
                        bbox_inches='tight',
                        dpi=300)
            pyplot.close(fig)




    # not used
    @staticmethod



    # adp.analysis
    @staticmethod
    def violinplot(
        filename,
        y_vals,
        x_labels,
        title,
        x_lab='x',
        y_lab='y',
        x_lim=None,
        y_lim=None,
        rotate_x_labels=True,
        hlines=[],
        vlines=[],
        ):
        """Generate standard violin plot"""

        fig, axis = pyplot.subplots(nrows=1, ncols=1)
        axis.set_title(title)
        axis.violinplot(y_vals, showmeans=True)
        axis.set_xticks(range(1,len(x_labels)+1))
        axis.set_xticklabels(x_labels)
        for v in hlines: axis.axhline(v)
        for v in vlines: axis.axvline(v)
        axis.set_xlabel(x_lab)
        axis.set_ylabel(y_lab)
        axis.set_xlim(x_lim)
        axis.set_ylim(y_lim)
        if rotate_x_labels:
            pyplot.setp(axis.get_xticklabels(), rotation=90)
        fig.tight_layout()
        fig.savefig(filename)#, dpi=300)
        pyplot.close(fig)

        return

    # not used
    @staticmethod
    def multi_hierarchy_bar(
        prefix,
        hierarchies,
        titles,
        shape=None,
        y_labs='Isotropic B ($\AA^2$)',
        v_line_hierarchy=None,
        rotate_x_labels=True,
        ):
        """Plot multiple bar plots for a series of hierarchies (plotted values are the average B-factors of each residue of the hierarchies)"""

        # Check titles same lengths
        titles = list(titles)
        assert len(hierarchies) == len(titles)
        # Check y-labs same lengths
        if isinstance(y_labs, str):
            y_labs = [y_labs]*len(hierarchies)
        assert len(hierarchies) == len(y_labs)
        # Get and check the shape of the subplots
        if shape is not None:
            nrow, ncol = shape
        else:
            nrow, ncol = (len(hierarchies), 1)
        assert len(hierarchies) == nrow*ncol

        # Extract one of the hierarchies as the "master"
        m_h = hierarchies[0]
        # Check all hierarchies are the same
        for h in hierarchies:
            assert m_h.is_similar_hierarchy(h)

        # Create a plot for each chain
        for chain_id in sorted(set([c.id for c in m_h.chains()])):

            # Create a selection for each chain and select it in each hierarchy
            sel = m_h.atom_selection_cache().selection('chain {}'.format(chain_id))
            sel_hs = [h.select(sel) for h in hierarchies]
            sel_mh = sel_hs[0]

            # Extract y-vals for each structure
            y_vals = []
            for i,h in enumerate(sel_hs):
                y_vals.append([numpy.mean(rg.atoms().extract_b()) for rg in h.residue_groups()])
            y_vals = numpy.array(y_vals)
            # Skip chains with no Bs
            if not y_vals.any():
                continue

            # Filename!
            filename = prefix + '-chain_{}.png'.format(chain_id)

            # Create x-values for each residue starting from 1
            x_vals = numpy.array(range(len(list(sel_mh.residue_groups()))))+1
            x_labels = ['']+[ShortLabeller.format(rg) for rg in sel_mh.residue_groups()]

            # Create the output figure
            fig, axes = pyplot.subplots(nrows=nrow, ncols=ncol)
            # Fix for 1 subplot
            if nrow==ncol==1: axes = numpy.array([axes])

            # Iterative though axes and plot
            handles = []
            for i, axis in enumerate(axes):
                # Set labels & title
                axis.set_ylabel(y_labs[i])
                axis.set_title(label=titles[i])
                # Plot the bar
                hdl = axis.bar(left=x_vals,
                               height=y_vals[i],
                               width=1.0,
                               align='center',
                               edgecolor='black')
                # Plot boundaries
                if v_line_hierarchy is not None:
                    v_lines = numpy.where(numpy.array([max(rg.atoms().extract_b()) for rg in v_line_hierarchy.select(sel).residue_groups()], dtype=bool))[0] + 1.5
                    for val in v_lines: axis.axvline(x=val, ls='dotted')

            # Make sure only the last axes have x-labels
            for a in axes[-ncol:]:
                a.set_xlabel('Residue')
                a.set_xticklabels([x_labels[int(i)] if (i<len(x_labels)) and (float(int(i))==i) else '' for i in axis.get_xticks()])

            # Format and save
            fig.tight_layout()
            fig.savefig(filename)#, dpi=300)
            pyplot.close(fig)

    # not used
    @classmethod
    def multi_scatter(cls, 
        filename,
        x=None, 
        x_vals=None, 
        y=None, 
        y_vals=None,
        legends=None,
        title='', 
        x_lab='x', 
        y_lab='y',
        x_lim=None, 
        y_lim=None,
        shape=None,
        hlines=[], 
        vlines=[],
        ):
        """Generate a scatter plot from data (or list of data)"""

        assert [(x is None),(x_vals is None)].count(False) == 1, 'must provide x OR x_vals'
        assert [(y is None),(y_vals is None)].count(False) == 1, 'must provide y OR y_vals'
        if x is not None: x_vals = [x]
        if y is not None: y_vals = [y]
        x_vals = [numpy.array(x) for x in x_vals]
        y_vals = [numpy.array(y) for y in y_vals]
        n_x = len(x_vals)
        n_y = len(y_vals)
        # Share X- or Y-axes?
        share_x = (n_x==1)
        share_y = (n_y==1)
        # Allow for differently sized lists if one is of unitary length
        #assert (n_x==1) or (n_y==1) or (n_x==n_y), 'number of x_vals not equal to number of y_vals'
        # Create cycles for simplicity
        x_vals = itertools.cycle(x_vals)
        y_vals = itertools.cycle(y_vals)
        # Create cycles of the x/y-labels
        x_lab = itertools.cycle([x_lab]) if isinstance(x_lab, str) else itertools.cycle(x_lab)
        y_lab = itertools.cycle([y_lab]) if isinstance(y_lab, str) else itertools.cycle(y_lab)

        # Create one plot if no shape specified
        if shape is None: shape = (1,1)

        n_line = max(n_x, n_y)
        n_plot = numpy.product(shape)

        if legends is not None:
            assert len(legends) == n_line

        # Colours of each of y_vals
        if n_plot == n_line:
            colours = ['b']*n_line
        elif (n_line % n_plot) == 0:
            colours = sorted(numpy.concatenate([pyplot.cm.rainbow(numpy.linspace(0,1,n_line//n_plot))]*n_plot).tolist())
        else:
            colours = pyplot.cm.rainbow(numpy.linspace(0,1,n_line))
        assert len(colours)==n_line

        # Create figures
        fig, axes = pyplot.subplots(nrows=shape[0], ncols=shape[1],
                                    sharex=share_x, sharey=share_y)
        # Create list if only one plot
        if shape == (1,1): axes = [axes]
        else: axes = numpy.array(axes).flatten()
        # Set the figure title
        fig.suptitle(title)
        # Draw horizontal/vertical lines (do first so they're at the bottom)
        for axis in axes:
            for v in hlines: axis.axhline(y=v, linewidth=2, zorder=1)
            for v in vlines: axis.axvline(x=v, linewidth=2, zorder=1)
        # Store plot objects
        plot_dicts = []
        axes_cycle = itertools.cycle(axes)
        for i in xrange(n_line):
            axis = axes_cycle.next()
            plt = axis.scatter(x=x_vals.next(),
                               y=y_vals.next(),
                               c=colours[i])
            if legends: plt.set_label(legends[i])
            plot_dicts.append(plt)
        # Format axes
        for axis in axes:
            # Axis labels
            axis.set_xlabel(x_lab.next())
            axis.set_ylabel(y_lab.next())
            # Make axis legends
            if legends:
                axis.legend()
        # Axis limits
        if x_lim: axis.set_xlim(x_lim)
        if y_lim: axis.set_ylim(y_lim)
        # X-axis rotations
        #if rotate_x_labels:
        #    pyplot.setp(axis.get_xticklabels(), rotation=45)
        fig.tight_layout()
        fig.savefig(filename)#, dpi=300)
        pyplot.close(fig)

        return

