import giant.logs as lg
logger = lg.getLogger(__name__)

import itertools, collections
import numpy

from libtbx import adopt_init_args, group_args

from giant.structure.formatting import ShortLabeller

try:
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.interactive(False)
    from matplotlib import pyplot
    pyplot.switch_backend('agg') # yes I know this done twice -- for safety!
    pyplot.interactive(0)
except Exception as e:
    logger(e)
    import matplotlib
    from matplotlib import pyplot

# def dendrogram(fname, link_mat, labels=None, ylab=None, xlab=None, ylim=None, annotate_y_min=0.25, num_nodes=20):
#     from matplotlib import pyplot
#     fig = pyplot.figure()
#     ax1 = pyplot.subplot(1,1,1)
#     dend = scipy.cluster.hierarchy.dendrogram(link_mat, p=num_nodes, truncate_mode='lastp', labels=labels)
#     # Change labels if requested
#     if xlab:   ax1.set_xlabel(xlab)
#     if ylab:   ax1.set_ylabel(ylab)
#     if ylim:   ax1.set_ylim(ylim)
#     # Make sure the labels are rotated
#     xlocs, xlabels = pyplot.xticks()
#     pyplot.setp(xlabels, rotation=90)
#     for i, d in zip(dend['icoord'], dend['dcoord']):
#         x = 0.5 * sum(i[1:3])
#         y = d[1]
#         if y < annotate_y_min: continue
#         pyplot.plot(x, y, 'ro')
#         pyplot.annotate("%.3g" % y, (x, y), xytext=(0, -8), textcoords='offset points', va='top', ha='center')
#     pyplot.tight_layout()
#     fig.savefig(fname)
#     return fig


class PiecewiseLinearFunction(object):

    def __init__(self, v_min, v_max, m, c):
        adopt_init_args(self, locals())

    def __call__(self, x):
        return max(self.v_min, min(self.v_max, self.m*x+self.c))


class PlotHelper(object):

    _ms_map = dict(
        # (min, max, m, c)
        narrow = PiecewiseLinearFunction(1.0, 2.0, -0.04, 2.8),        # (1, 3, -0.1, 5)
        chunky = PiecewiseLinearFunction(2.0, 4.0, -0.08, 5.6),        # (3, 5, -0.1, 7)
        )
    _lw_map = dict(
        narrow = PiecewiseLinearFunction(0.5, 1.0, -0.02, 1.4),     # (1,1,0,0)
        chunky = PiecewiseLinearFunction(1.0, 2.0, -0.04, 2.8),     # (2,2,0,0)
        )
    _labelsize_map = PiecewiseLinearFunction(6, 10, -0.15, 14)

    default = group_args(
        marker = 'D',
        )

    def __init__(self,
        colour_map_name = 'rainbow',
        plot_style = 'ggplot',
        font_family = 'monospace',
        font_name = None,
        dpi = 300,
        ):

        adopt_init_args(self, locals())

        if (plot_style is not None) and (plot_style not in ['xkcd']):
            try:
                logger('Setting plot_style to "{}"'.format(plot_style))
                pyplot.style.use(plot_style)
            except Exception as e:
                logger.warning('Failed to set plot style to {}.\n\t{}'.format(plot_style, str(e)))

        if (plot_style == 'xkcd'):
            try:
                pyplot.xkcd()
                logger('Theme set to "xkcd".')
            except:
                logger.warning('Failed to set plot style "xkcd".')

        if (font_family is not None):
            try:
                logger('Setting font_family to "{}"'.format(font_family))
                pyplot.rc('font', family=font_family)
            except Exception as e:
                logger.warning('Failed to set font family to {}.\n\t{}'.format(font_family, str(e)))

        if (font_name is not None):
            try:
                assert font_family is not None, 'Cannot set font: must provide a font family in order to set font.'
                font_family_str = 'font.'+str(font_family)
                assert font_family_str in list(pyplot.rcParams.keys()), 'Cannot set font: invalid font family provided "{}".'.format(font_family)
                family_fonts = pyplot.rcParams[font_family_str]
                if (font_name not in family_fonts):
                    logger.warning('WARNING: font "{}" does not exist in font family "{}". Setting the font may not work... (valid options: {})'.format(font_name, font_family, ', '.join(family_fonts)))
                logger('Setting font_name to "{}"'.format(font_name))
                pyplot.rcParams[font_family_str].insert(0, font_name)
            except Exception as e:
                logger.warning('Failed to set font to {}.\n\t{}'.format(font_name, str(e)))

    def get_colour_map(self):
        return matplotlib.cm.get_cmap(self.colour_map_name)

    def get_colours(self, n):
        cm = self.get_colour_map()
        return cm(numpy.linspace(0., 1., n))

    def ms(self, n_x_values, ms_type='narrow'):
        return self._ms_map[ms_type](n_x_values)

    def lw(self, n_x_values, lw_type='narrow'):
        return self._lw_map[lw_type](n_x_values)

    def labelsize(self, n_x_values):
        return self._labelsize_map(n_x_values)

    def initialise_figure(self,
        title,
        x_label,
        y_label,
        axis = None,
        ):

        if axis is None:
            fig, axis = pyplot.subplots(nrows=1, ncols=1)
        else:
            fig = axis.get_figure()

        self.set_axis_labels(
            axis = axis,
            title = title,
            x_label = x_label,
            y_label = y_label,
            )
        return fig, axis

    def set_axis_labels(self,
        axis,
        title,
        x_label,
        y_label,
        ):
        # Set the font size to the "maximum" font size
        fd = dict(fontsize=self.labelsize(0))
        if title is not None:
            axis.set_title(title, fontdict=fd)
        if x_label is not None:
            axis.set_xlabel(x_label, fontdict=fd)
        if y_label is not None:
            axis.set_ylabel(y_label, fontdict=fd)

    def make_x_ticks(self,
        axis,
        x_ticks,
        x_tick_labels = None,
        n_labels = 10,
        ):
        n_total = len(x_ticks)
        if (x_tick_labels is None):
            x_tick_labels = list(map(str, x_ticks))
        i_x_ticks = numpy.arange(0, n_total, int(max(1.0, numpy.floor(float(n_total)/float(n_labels)))))
        axis.set_xticks(
            [x_ticks[i] for i in i_x_ticks],
            )
        axis.set_xticklabels(
            [x_tick_labels[i] for i in i_x_ticks],
            fontdict = dict(fontsize=self.labelsize(len(i_x_ticks))),
            )

    def rotate_x_tick_labels(self, axis, angle=90):
        """Convenience function that doesn't require importing pyplot"""
        pyplot.setp(axis.get_xticklabels(), rotation=angle)

    def hide_x_labels(self, axis):
        pyplot.setp(axis.get_xticklabels(), visible=False) # doesn't always work...

    def write_and_close_fig(self, fig, filename, **kw_args):
        assert fig is not None
        fig.tight_layout()
        fig.savefig(filename, bbox_inches='tight', dpi=self.dpi, **kw_args)
        pyplot.close(fig)

    def bin_x_values(self,
        data,
        n_bins=10,
        ):
        """Generate bins for data and return binned values and indices of binning"""

        data = numpy.array(data)
        bins = numpy.linspace(data.min(), data.max(), n_bins+1)
        # To avoid extra bin for the max value
        bins[-1] += 0.001*abs(bins[-1])
        # Assign data to bins
        indices = numpy.digitize(data, bins)
        # Generate bin labels
        bin_labels = ['{:.2f} - {:.2f}'.format(bins[i],bins[i+1]) for i in range(n_bins)]
        return bins, indices, bin_labels

    def resolve_value_arrays(self,
        x_vals,
        x_vals_array,
        y_vals,
        y_vals_array,
        ):
        assert [x_vals is None, x_vals_array is None].count(True) == 1
        assert [y_vals is None, y_vals_array is None].count(True) == 1

        if (x_vals is not None) and (y_vals is not None):
            x_vals_array = [x_vals]
            y_vals_array = [y_vals]
        elif (x_vals is None) and (y_vals is None):
            assert len(x_vals_array) == len(y_vals_array)
        elif (x_vals is not None):
            x_vals_array = [x_vals] * len(y_vals_array)
        else:
            y_vals_array = [y_vals] * len(x_vals_array)

        assert len(x_vals_array) == len(y_vals_array)
        n = len(x_vals_array)

        return (n, x_vals_array, y_vals_array)


class PandemicAdpPlotter(object):
    """
    Generic plotting class. Needs to be initialised for some functions but not for all.
    """

    helper = None

    function_hash = {
        'min' : (numpy.min, 'k:'),
        'mean': (numpy.mean, 'k-'),
        'max' : (numpy.max, 'k:'),
    }

    def __init__(self,
        n_levels = None,
        helper = None,
        ):

        adopt_init_args(self, locals(), exclude=('helper',))

        # Use provided if given
        if (helper is not None):
            self.helper = helper
        # If class variable is not initialised then initialise new instance
        elif (self.helper is None):
            self.helper = PlotHelper()

    def get_level_colours_arbitrary(self, indices):
        cm = self.helper.get_colour_map()
        return cm(numpy.array(indices)/float(self.n_levels-1))

    def get_level_colours(self):
        cm = self.helper.get_colour_map()
        return cm(numpy.linspace(0., 1., self.n_levels))

    def multi_hierarchy_plot_by_residue(self,
        hierarchies,
        plot_function,
        plot_kw_args,
        prefix,
        residue_values_function = None,
        y_array_values_function = None,
        ):
        """Plot a graph binned by each residue of a structure, where the input values are the B-factors of the atoms. Produces one image for each chain."""

        if residue_values_function is None:
            residue_values_function = self.array_passthrough

        if y_array_values_function is None:
            y_array_values_function = self.array_passthrough

        output_files = collections.OrderedDict()

        m_h = hierarchies[0]

        # Create a plot for each chain
        for chain_id in sorted(set([c.id for c in m_h.chains()])):

            # Create a selection for each chain and select it in each hierarchy
            sel = m_h.atom_selection_cache().selection('chain {}'.format(chain_id))
            hierarchies_sel = [h.select(sel) for h in hierarchies]
            m_h_sel = m_h.select(sel)

            # Filename!
            filename = prefix + '-chain_{}.png'.format(chain_id)
            output_files[chain_id] = filename

            # Create x-values for each residue starting from 1
            x_tick_labels = [ShortLabeller.format(rg) for rg in m_h_sel.residue_groups()]
            x_vals = 1.0 + numpy.arange(len(x_tick_labels))

            # Mean values to be plotted as line
            y_vals_array = []
            for h in hierarchies_sel:
                y_vals = [residue_values_function(list(rg.atoms().extract_b())) for rg in h.residue_groups()]
                y_vals_array.append(y_vals)

            # Post-process the extracted values (for instance, combine all hierarchies into one using array_concatentate)
            y_vals_array = y_array_values_function(y_vals_array)
            n = max(list(map(len,y_vals_array)))

            #
            kw_args = {}
            kw_args.update(plot_kw_args)

            plot_function(
                x_vals = x_vals,
                y_vals_array = y_vals_array,
                x_ticks = x_vals,
                x_tick_labels = x_tick_labels,
                filename = filename,
                **plot_kw_args
                )

        return output_files

    @staticmethod
    def array_passthrough(values):
        return values

    @staticmethod
    def array_concatenate(values):
        new_values = []
        for old_values in zip(*values):
            l = numpy.concatenate([numpy.array(v).flatten() for v in old_values])
            new_values.append(l)
        return [new_values]

    def scatter(self,
        axis = None,
        x_vals = None,
        x_vals_array = None,
        y_vals = None,
        y_vals_array = None,
        title = '',
        x_label = '',
        y_label = '',
        x_ticks = None,
        x_tick_labels = None,
        x_lim = None,
        y_lim = None,
        alphas = None,
        legends = None,
        filename = None,
        **plot_kw_args
        ):

        fig, axis = self.helper.initialise_figure(
            title = title,
            x_label = x_label,
            y_label = y_label,
            )

        n, x_vals_array, y_vals_array = self.helper.resolve_value_arrays(
            x_vals = x_vals,
            x_vals_array = x_vals_array,
            y_vals = y_vals,
            y_vals_array = y_vals_array,
            )

        colours = self.helper.get_colours(n)

        if alphas is None:
            alphas = [1.] * n
        elif isinstance(alphas, float):
            alphas = [alphas] * n
        else:
            assert len(alphas) == n

        handles = []
        for i, (x, y) in enumerate(zip(x_vals_array, y_vals_array)):
            kw_args = {
                'marker' : self.helper.default.marker,
                'color' : colours[i],
                'alpha' : alphas[i],
                }
            kw_args.update(plot_kw_args)
            hdl = axis.scatter(
                x = x,
                y = y,
                **kw_args
                )
            handles.append(hdl)

        if (x_ticks is not None):
            self.helper.make_x_ticks(
                axis = axis,
                x_ticks = x_ticks,
                x_tick_labels = x_tick_labels,
                n_labels = 20,
                )

        if y_lim is not None:
            axis.set_ylim(y_lim)

        if legends is not None:
            assert len(legends) == n
            assert len(handles) == n
            axis.legend(
                handles = handles,
                labels = legends,
                fontsize = self.helper.labelsize(0),
                )

        if filename is not None:
            self.helper.write_and_close_fig(fig=fig, filename=filename)

        return fig, axis

    def lineplot(self,
        axis = None,
        x_vals = None,
        x_vals_array = None,
        y_vals = None,
        y_vals_array = None,
        title = '',
        x_label = '',
        y_label = '',
        x_ticks = None,
        x_tick_labels = None,
        y_lim = None,
        rotate_x_labels = True,
        alphas = None,
        legends = None,
        filename = None,
        background_line_type = None,
        legend_kw_args = dict(bbox_to_anchor=(1.0, -0.15), loc=1, borderaxespad=0.),
        **plot_kw_args
        ):

        fig, axis = self.helper.initialise_figure(
            title = title,
            x_label = x_label,
            y_label = y_label,
            axis = axis,
            )

        n, x_vals_array, y_vals_array = self.helper.resolve_value_arrays(
            x_vals = x_vals,
            x_vals_array = x_vals_array,
            y_vals = y_vals,
            y_vals_array = y_vals_array,
            )

        colours = self.helper.get_colours(n)

        if alphas is None:
            alphas = [1.] * n
        elif isinstance(alphas, float):
            alphas = [alphas] * n
        else:
            assert len(alphas) == n

        handles = []
        for i, (x, y) in enumerate(zip(x_vals_array, y_vals_array)):
            nx = len(x)
            kw_args = {
                'marker' : self.helper.default.marker,
                'color' : colours[i],
                'alpha' : alphas[i],
                'linewidth' : self.helper.lw(nx),
                'markersize' : self.helper.ms(nx),
                }
            kw_args.update(plot_kw_args) # allow overrides of defaults

            # Plot a black line underneath the coloured lines
            if background_line_type is not None:
                bg_kw_args = dict(**kw_args)
                # override any custom provided args
                bg_kw_args['linewidth'] = self.helper.lw(nx, background_line_type)
                bg_kw_args['markersize'] = self.helper.ms(nx, background_line_type)
                bg_kw_args['color'] = 'k'
                hdl = axis.plot(x, y, **bg_kw_args)

            hdl = axis.plot(x, y, **kw_args)
            handles.extend(hdl)

        if (x_ticks is not None):
            self.helper.make_x_ticks(
                axis = axis,
                x_ticks = x_ticks,
                x_tick_labels = x_tick_labels,
                n_labels = 20,
                )

        if y_lim is not None:
            axis.set_ylim(y_lim)

        artists = []

        if legends is not None:
            assert len(legends) == n
            assert len(handles) == n
            lgd = axis.legend(
                handles = handles,
                labels = legends,
                fontsize = self.helper.labelsize(0),
                **legend_kw_args
                )
            artists.append(lgd)

        if rotate_x_labels:
            pyplot.setp(axis.get_xticklabels(), rotation=90)

        if filename is not None:
            self.helper.write_and_close_fig(
                fig = fig,
                filename = filename,
                bbox_extra_artists = artists,
                )

        return fig, axis

    def boxplot(self,
        x_vals = None,
        x_vals_array = None,
        y_vals = None,
        y_vals_array = None,
        title = '',
        x_label = '',
        y_label = '',
        x_ticks = None,
        x_tick_labels = None,
        y_lim = None,
        hlines = [],
        vlines = [],
        rotate_x_labels = True,
        legends = None,
        filename = None,
        **plot_kw_args
        ):
        """Generate standard boxplot"""

        assert x_vals_array is None

        fig, axis = self.helper.initialise_figure(
            title = title,
            x_label = x_label,
            y_label = y_label,
            )

        if (x_vals is None):
            if (y_vals is not None):
                x_vals = 1 + numpy.arange(len(y_vals))
            else:
                x_vals = 1 + numpy.arange(len(y_vals_array[0]))

        n, x_vals_array, y_vals_array = self.helper.resolve_value_arrays(
            x_vals = x_vals,
            x_vals_array = x_vals_array,
            y_vals = y_vals,
            y_vals_array = y_vals_array,
            )

        if (x_ticks is None):
            x_ticks = x_vals

        if (x_tick_labels is None):
            x_tick_labels = list(map(str, x_ticks))

        for i, (x, y) in enumerate(zip(x_vals_array, y_vals_array)):
            kw_args = {
                'meanprops' : {
                    'marker' : self.helper.default.marker,
                    'markersize' : self.helper.ms(len(x)),
                    },
                }
            kw_args.update(plot_kw_args)
            axis.boxplot(
                x = y,
                positions = x,
                labels = x_tick_labels,
                showmeans = True,
                **kw_args
                )

        if (x_ticks is not None):
            self.helper.make_x_ticks(
                axis = axis,
                x_ticks = x_ticks,
                x_tick_labels = x_tick_labels,
                n_labels = 20,
                )

        if (y_lim is not None):
            axis.set_ylim(y_lim)

        for l in hlines: axis.axhline(l)
        for l in vlines: axis.axvline(l)

        if rotate_x_labels:
            pyplot.setp(axis.get_xticklabels(), rotation=90)

        if filename is not None:
            self.helper.write_and_close_fig(fig=fig, filename=filename)

        return fig, axis

    # adp.process_output - REFACTOR
    def binned_boxplot(self,
        filename,
        x,
        y=None,
        y_vals=None,
        legends=None,
        title='',
        x_lab='x',
        y_lab='y',
        rotate_x_labels=True,
        plot_scatter=True,
        max_bins=10,
        min_bin_width=None,
        hlines=[],
        vlines=[],
        plot_type='violinplot',
        ):
        """Generate a binned boxplot from data (or array of data)"""

        assert [(y is None),(y_vals is None)].count(False) == 1, 'must provide y OR y_vals'
        if y is not None: y_vals = [y]
        y_vals = [numpy.array(y) for y in y_vals]
        n_y = len(y_vals)

        if legends is not None:
            assert len(legends) == n_y

        assert isinstance(max_bins, int) and max_bins>0

        if min_bin_width is not None:
            assert min_bin_width > 0
            min_x = numpy.min(x)
            max_x = numpy.max(x)
            n_bins = int(min(max_bins, 1+((max_x-min_x)//min_bin_width)))
        else:
            n_bins = max_bins
        assert n_bins >= 1

        # Stupid list of things to not cut off of the graph
        extra_artists = []

        # Generate binning for the x axis
        bins, indices, bin_labels = self.helper.bin_x_values(data=x, n_bins=n_bins)
        # Sort the y_vals into the bins
        binned_y = [[y[indices==i] for i in range(1, n_bins+1)] for y in y_vals]
        # Width of the boxplot bars
        bar_width = 2./float(1+3*n_y) # half-bar padding between bars
        # Colours of each of y_vals
        colours = pyplot.cm.rainbow(numpy.linspace(0,1,n_y))

        # Create figures
        fig, axis = pyplot.subplots(nrows=1, ncols=1)
        axis.set_title(title)
        # Draw horizontal/vertical lines (do first so they're at the bottom)
        for v in hlines: axis.axhline(y=v, linewidth=1, zorder=1)
        for v in vlines: axis.axvline(x=v, linewidth=1, zorder=1)
        # Store plot objects
        plot_dicts = []
        for i_y, y in enumerate(binned_y):
            # Offset for this bar set relative to (n-1)
            x_offset = 0.5 + (1+1.5*i_y)*bar_width
            positions = numpy.arange(n_bins) + x_offset
            # Filter on values that are actually present
            y_idx = [i for i in range(len(y)) if len(y[i])>0]
            # Plot
            if plot_type == 'boxplot':
                plt = axis.boxplot(
                        y,
                        positions=positions,
                        widths=bar_width,
                        showmeans=False,
                        patch_artist=True)
            else:
                plt = axis.violinplot(
                        [y[i] for i in y_idx],
                        positions=[positions[i] for i in y_idx],
                        widths=bar_width,
                        showmeans=True,
                        showextrema=True)
            plot_dicts.append(plt)
            # Color the boxplots
            c = colours[i_y]
            for obj in ['boxes','bodies']:
                for patch in plt.get(obj, []):
                    patch.set_facecolor(c)
                    patch.set_edgecolor('k')
            for obj in ['means','medians']:
                for line in plt.get(obj, []):
                    line.set_facecolor(c)
                    line.set_edgecolor('k')
            for obj in ['cbars','cmeans','cmedians','cmins','cmaxes']:
                lc = plt.get(obj, None)
                if lc is not None:
                    lc.set_facecolor(c)
                    lc.set_edgecolor('k')
            # Plot a superposed scatter plot
            if plot_scatter is True:
                for x, ys in zip(positions, y):
                    if len(ys) == 0: continue
                    xs = [x]*len(ys)
                    jitter = numpy.random.randn(len(ys)) * 0.5 * bar_width / 3.0 # Normalise to be approx width of bar
                    jitter = jitter - jitter.mean()
                    axis.scatter(
                        x = xs+jitter,
                        y = ys,
                        s = self.helper.ms(len(positions)*len(xs)),
                        facecolor = c,
                        edgecolor = 'k',
                        zorder = 10,
                        )

        # Make labels
        bin_centres = numpy.arange(n_bins)+1
        axis.set_xticks(bin_centres)
        axis.set_xticklabels(bin_labels)
        # Axis labels
        axis.set_xlabel(x_lab)
        axis.set_ylabel(y_lab)
        # Axis limits
        axis.set_xlim((0.25, n_bins+0.75))
        # Plot v_lines
        if n_y > 1:
            for i in range(n_bins+1):
                axis.axvline(i+0.5, c='grey', ls="dashed", lw=1.0)#ls="solid", lw=0.5)
        # X-axis rotations
        if rotate_x_labels:
            pyplot.setp(axis.get_xticklabels(), rotation=45)
        # Make legend
        if legends is not None:
            # Plot additional text on the graph
            y_min, y_max = axis.get_ylim()
            fig_scale = (axis.transData.transform((1.0,0.)) - axis.transData.transform((0.0,0.)))[0]
            fig_width = 30.0 * n_y # final width of the fonts in the scale of the image
            x_width = fig_width / fig_scale # Width on the x-axis scale
            t_width = x_width / float(n_y)
            fontsize = 10
            text_artists = []
            max_l_length = max(map(len,legends))
            for i_y, l in enumerate(legends):
                y_pos = (0.5*y_min+0.5*y_max) #(0.05*(y_max-y_min) if (y_max>0.0>y_min) else (0.5*y_min+0.5*y_max))
                t = axis.text(
                        x=(n_bins+0.5)+(0.5+i_y)*t_width, y=y_pos, s=l.center(max_l_length),
                        bbox=dict(boxstyle='round', facecolor=colours[i_y], edgecolor='k', linewidth=0.5, alpha=0.75, pad=0.3),
                        fontsize=fontsize, rotation=90, ha='center', va='center', zorder=1)
                text_artists.append(t)
            axis.set_xlim(right=n_bins+0.5+x_width)

        self.helper.write_and_close_fig(
            fig = fig,
            filename = filename,
            bbox_extra_artists = extra_artists+text_artists,
            )

        return

    # echt.process_output - REFACTOR
    def multi_histogram(self,
        filename,
        x_vals,
        titles,
        x_labs,
        rotate_x_labels=True,
        shape=None,
        n_bins=30,
        x_lim=None,
        ):
        """Generate standard histogram"""

        if shape is not None:
            nrow, ncol = shape
        else:
            nrow, ncol = (len(x_vals), 1)
        assert nrow*ncol == len(x_vals)

        fig, axes = pyplot.subplots(nrows=nrow, ncols=ncol, sharey=False)

        # Fix for 1 subplot
        if nrow==ncol==1: axes = numpy.array([axes])

        for i, axis in enumerate(axes.flatten()):

            # Select x-axis limits
            if x_lim is not None:
                assert len(x_lim) == 2
                i_x_lim = [min(numpy.min(x_vals[i]), x_lim[0]) if (x_lim[0] is not None) else numpy.min(x_vals[i]),
                           max(numpy.max(x_vals[i]), x_lim[1]) if (x_lim[1] is not None) else numpy.max(x_vals[i])]
            else:
                i_x_lim = [numpy.min(x_vals[i]), numpy.max(x_vals[i])]
            # Make sure does not have zero width
            if (i_x_lim[1] - i_x_lim[0]) < 1e-6:
                i_x_lim[0] -= 0.01*abs(i_x_lim[0])
                i_x_lim[1] += 0.01*abs(i_x_lim[1])

            axis.set_title(titles[i])
            axis.hist(x=x_vals[i], bins=n_bins, range=i_x_lim)
            axis.set_xlabel(x_labs[0])
            axis.set_ylabel('Count')
            axis.set_xlim(i_x_lim)
            if rotate_x_labels:
                pyplot.setp(axis.get_xticklabels(), rotation=90)
        fig.tight_layout()
        fig.savefig(filename)#, dpi=300)
        pyplot.close(fig)

        return

    # adp.hierarchy
    def level_plots(self,
        filename,
        hierarchies,
        labels,
        title,
        rotate_x_labels=True,
        ):
        """Plot a schematic representation of the hierarchical partitioning"""

        colours = self.get_level_colours()

        fig, axis = self.helper.initialise_figure(
            title = title,
            x_label = 'Atom',
            y_label = '',
            )
        axis.set_facecolor('w')

        assert len(hierarchies) <= (self.n_levels)

        for i_h, h in enumerate(hierarchies):
            # Extract B-factors and atom labels
            b_vals = numpy.array(h.atoms().extract_b())
            n_atoms = len(b_vals)
            # Extract colour for this level
            col = colours[i_h]
            # Plot a bar around the atoms (below other bars)
            axis.broken_barh([(0.25, n_atoms+0.50)], (i_h+0.5, 1.0), edgecolor=None, facecolor='w')
            # Iterate through b values and draw boxes for each
            for b in numpy.unique(b_vals):
                # Skip zero or negative values
                if b < 0: continue
                # Indices of values in this group
                idx = numpy.where((b_vals==b))[0]
                # Cluster the indices to create bars
                if len(idx)==1:
                    cluster = [1]
                else:
                    import scipy.cluster
                    cluster = scipy.cluster.hierarchy.fclusterdata(X=idx.reshape((len(idx),1)), t=1.1, criterion='distance', metric='euclidean', method='single')
                # Iterate through clusters and draw
                plot_vals = []
                for g in numpy.unique(cluster):
                    c_idxs = idx[g==cluster]
                    min_i = numpy.min(c_idxs)
                    max_i = numpy.max(c_idxs)
                    plot_vals.append((min_i+0.5, max_i-min_i+1.00))
                # Plot
                axis.broken_barh(plot_vals, (i_h+0.5, 1.0), edgecolor='k', facecolor=col)

        # Set axes, etc.
        a_labels = numpy.array([ShortLabeller.format(a.parent().parent()) for a in hierarchies[0].atoms()])

        self.helper.make_x_ticks(
            axis = axis,
            x_ticks = 1.0 + numpy.arange(len(a_labels)),
            x_tick_labels = a_labels,
            n_labels = 30,
            )

        axis.set_yticks(list(range(1, len(hierarchies)+1)))
        axis.set_yticklabels(labels)
        axis.set_xlim((-0.01*len(a_labels), 1.01*(len(a_labels)+1)))
        pyplot.setp(axis.get_xticklabels(), rotation=90)
        for l in axis.get_yticklabels():
            l.set_multialignment('right')
        # Format and save
        self.helper.write_and_close_fig(fig=fig, filename=filename)

    # echt.process_output
    def stacked_bar(self,
        prefix,
        hierarchies,
        legends,
        reference_hierarchy=None,
        reference_legend=None,
        reference_functions=['mean'],
        title=None,
        y_lab='Isotropic B ($\\AA^2$)',
        y_lim=None,
        v_line_hierarchy=None,
        rotate_x_labels=True,
        reverse_legend_order=False,
        colour_space=(0,1),
        colour_indices=None,
        colours=None,
        legend_kw_args = dict(bbox_to_anchor=(1.0, 0.0), loc=1, borderaxespad=0.),
        ):
        """Plot stacked bar plots for a series of hierarchies (plotted values are the average B-factors of each residue of the hierarchies)"""

        legends = list(legends)
        assert len(hierarchies) == len(legends)

        m_h = hierarchies[0]

        if colour_indices is not None:
            assert len(hierarchies) == len(colour_indices)
            colours = self.get_level_colours_arbitrary(colour_indices)
        elif colours is not None:
            assert len(hierarchies) == len(colours)
        else:
            cm = self.helper.get_colour_map()
            colours = cm(numpy.linspace(colour_space[0],colour_space[1],len(hierarchies)))

        # Check all hierarchies are the same
        for h in hierarchies:
            assert m_h.is_similar_hierarchy(h)
        if reference_hierarchy is not None:
            assert m_h.is_similar_hierarchy(reference_hierarchy)

        output_files = collections.OrderedDict()

        # Create a plot for each chain
        for chain_id in sorted(set([c.id for c in m_h.chains()])):

            # Create a selection for each chain and select it in each hierarchy
            sel = m_h.atom_selection_cache().selection('chain {}'.format(chain_id))
            sel_hs = [h.select(sel) for h in hierarchies]
            sel_mh = sel_hs[0]

            # Filename!
            filename = prefix + '-chain_{}.png'.format(chain_id)
            output_files[chain_id] = filename

            # Create x-values for each residue starting from 1
            x_vals = numpy.array(list(range(len(list(sel_mh.residue_groups())))))+1
            x_labels = ['']+[ShortLabeller.format(rg) for rg in sel_mh.residue_groups()]
            # Cumulative y-values (for bottoms of bars)
            cuml_y = None

            # Create colours + hatches
            #hatchs = itertools.cycle(['//', 'x', '\\'])
            hatchs = itertools.cycle([None])

            # Create the output figure
            fig, axis = pyplot.subplots(nrows=1, ncols=1)
            if title is not None: axis.set_title(label=str(title))
            axis.set_xlabel('Residue')
            axis.set_ylabel(y_lab)
            if y_lim:
                axis.set_ylim(y_lim)

            # Iterative though hierarchies
            handles = []
            for i_h, h in enumerate(sel_hs):
                # Extract b-factors from this hierarchy
                y_vals = numpy.array([numpy.mean(rg.atoms().extract_b()) for rg in h.residue_groups()])
                # Initialise cumulative object
                if cuml_y is None:
                    cuml_y = numpy.zeros_like(y_vals)
                # Plot the bar
                hdl = axis.bar(
                        x=x_vals,
                        height=y_vals,
                        bottom=cuml_y,
                        width=1.0,
                        align='center',
                        color=colours[i_h],
                        label=legends[i_h],
                        hatch=next(hatchs),
                        linewidth=0,
                        edgecolor='black',
                        zorder=5,
                        )
                handles.append(hdl)

                # Append to cumulative y
                cuml_y += y_vals

            # Plot lines at top of graph for prettiness/as reference hierarchy
            if (reference_hierarchy is not None) or (len(sel_hs)==1):
                if reference_hierarchy is not None:
                    line_h_sel = reference_hierarchy.select(sel)
                else:
                    line_h_sel = sel_hs[0]
                # Plot one line for each function
                for ref_func in reference_functions:
                    # Extract the function from the hash
                    func, ls = self.function_hash[ref_func]
                    # Extract y-values
                    y_vals = numpy.array([func(rg.atoms().extract_b()) for rg in line_h_sel.residue_groups()])
                    # Duplicate the x-vals and y-vals and shift x-vals apart by half a unit
                    x_vals_dup = numpy.concatenate([(x_vals-0.5),(x_vals+0.5)]).reshape((2,x_vals.size)).T.reshape(2*x_vals.size).tolist()
                    y_vals_dup = numpy.concatenate([y_vals,y_vals]).reshape((2,y_vals.size)).T.reshape(2*y_vals.size).tolist()
                    # Add another point at the beginning and end so starts and ends on the baseline
                    x_vals_dup = numpy.concatenate([[x_vals_dup[0]], x_vals_dup, [x_vals_dup[-1]]])
                    y_vals_dup = numpy.concatenate([[0.0], y_vals_dup, [0.0]])
                    # Legend for line
                    ref_leg_label = reference_legend
                    hdl = axis.plot(x_vals_dup, y_vals_dup, ls, label=ref_leg_label, lw=0.5, zorder=5)
                    if reference_hierarchy is not None:
                        handles.extend(hdl)

            # Legends (reverse the plot legends!)
            if reverse_legend_order is True:
                handles.reverse()

            # No lines if cuml_y is None
            if cuml_y is None:
                continue

            # Plot boundaries
            if v_line_hierarchy is not None:
                v_lines = numpy.where(numpy.array([max(rg.atoms().extract_b()) for rg in v_line_hierarchy.select(sel).residue_groups()], dtype=bool))[0] + 1.5
                for val in v_lines:
                    axis.axvline(x=val, c='grey', ls='solid', label='boundaries', lw=0.5, zorder=1)
                # Add a line at zero
                h = axis.axvline(x=0.5, c='grey', ls='solid', label='boundaries', lw=0.5, zorder=1)
                handles.append(h)

            # Plot legend
            lgd = fig.legend(
                handles = handles,
                fontsize = self.helper.labelsize(0),
                bbox_transform=fig.transFigure,
                **legend_kw_args
                )

            # Axis ticks & labels
            x_ticks = numpy.arange(1, len(x_labels)+1, int(max(1.0, numpy.floor(float(len(x_labels))/20.0))))
            #x_ticks = numpy.unique(numpy.floor(numpy.linspace(1,len(x_labels)+1,20)))
            axis.set_xticks(x_ticks)
            axis.set_xticklabels([x_labels[int(i)] if (i<len(x_labels)) and (float(int(i))==i) else '' for i in axis.get_xticks()])
            # Rotate axis labels
            if rotate_x_labels: pyplot.setp(axis.get_xticklabels(), rotation=90)

            # Format and save
            self.helper.write_and_close_fig(
                fig = fig,
                filename = filename,
                bbox_extra_artists=[lgd],
                )

        return output_files

    @staticmethod
    def failure_graph(filename, exception=None, title=None, **args):
        fig = pyplot.figure()
        if title is None: title=filename
        if exception is not None: title += '\n---\n'+str(exception) + '\n---'
        pyplot.title('Failed to make {}'.format( title ))
        fig.savefig(filename)
        pyplot.close(fig)
