from matplotlib import pyplot
pyplot.style.use('ggplot')

import numpy

class Radar(object):

    title_fontsize = 12
    label_fontsize = 10
    legend_fontsize = 12

    fancy_labels = True

    # Inside and Outside offsets
    _padding_inner = 0.4
    _padding_outer = 0.3

    def __init__(self, titles, rect=None, fig=None):

        # Make a rectangle
        if rect is None:
            rect = [0.05, 0.05, 0.95, 0.95]
        # Number of plots
        self.n = len(titles)
        # Number of data points
        self.d = 0

        # Control flags
        self._plotted    = False

        # Axis limits + ranges + ticks
        self._limits = []
        self._ranges = []
        self._ticks  = []
        # Scaling to 0-1 (Scale all to this scale)
        self._scales = []
        self._offset = []
        # Data Points
        self._data = []
        self._args = []
        self._kwgs = []
        # Invert Axes?
        self._invert = [0]*self.n

        # Create a figure if none given
        if (fig is None):
            fig = pyplot.figure()
        self.fig = fig

        # Store titles
        self.titles = titles
        # Create and store angles and axes
        self.angles = [a if (a <= 360.) else (a - 360.) for a in numpy.arange(90, 90+360, 360.0/self.n)]

        # Create main plotting axis (all others are only used for axis labels)
        self.ax = self.fig.add_axes(
            rect,
            projection = "polar",
            label = "main-axis",
        )
        # Create individual axes for each spoke
        self.axes = [self.fig.add_axes(rect, projection="polar", label="axes%d" % i) for i in range(self.n)]

        # Background colour
        self.ax.set_facecolor('lightgrey')

        # Add grids and labels
        self.ax.set_thetagrids(
            self.angles,
            labels = ['']*self.n,
        )

        # Plot axis labels separately
        for i, (t, ang) in enumerate(zip(titles, numpy.deg2rad(numpy.r_[self.angles, self.angles[0]]))):
            ha, va = self._ha_va(i)
            self.ax.text(
                x = ang,
                y = self._padding_inner + 1.0 + self._padding_outer,
                s = t,
                horizontalalignment = ha,
                verticalalignment = va,
                fontsize = self.title_fontsize,
                weight = 'bold',
            )

        # Hide unused graph backgrounds
        self.ax.fill_between(
            numpy.deg2rad(numpy.r_[self.angles, self.angles[0]]),
            0, self._padding_inner,
            lw = 0,
            facecolor = 'white',
        )
        self.ax.fill_between(
            numpy.deg2rad(numpy.r_[self.angles, self.angles[0]]),
            self._padding_inner + 1., 10,
            lw = 0,
            facecolor = 'white',
        )

        self.ax.xaxis.grid(
            color = 'black',
            linestyle = '-',
            linewidth = 1,
        )
        self.ax.yaxis.grid(False)
        self.ax.yaxis.set_visible(False)

        # Remove all things from other axes
        for ax in self.axes:
            ax.patch.set_visible(False)
            ax.xaxis.set_visible(False)
            ax.grid(False)

    def _normalise(self, vals):
        """Normalise values to the scaled axis"""
        assert len(vals) == self.n
        ret_vals = []
        for vls,s,o,i in zip(vals, self._scales, self._offset, self._invert):
            ret_vals_sub = []
            for v in vls:
                if (v is None) or numpy.isnan(v):
                    ret_vals_sub.append(self._padding_inner+0.5)
                else:
                    ret_vals_sub.append(self._padding_inner+(1.0-i)*(v*s+o)+i*(1.0-(v*s+o)))
            ret_vals.append(ret_vals_sub)
        return ret_vals

    def _default_limits(self):
        """Set the limits to the range of the data"""
        limits = []
        for d in numpy.array(self._data).T:
            if (d==None).any():  limits.append((0.49, 0.51))
            elif len(d)==1:      limits.append((d[0]-0.01, d[0]+0.01))
            elif min(d)==max(d): limits.append((d[0]-0.01, d[0]+0.01))
            else:                limits.append((min(d), max(d)))
        self.set_limits(limits)

    def _calculate_scales(self):
        """Calculate the scales to normalise the axis limits to {0,1}"""
        # Calculate ranges for each axis
        self._ranges = [1.0*(max-min) for (min, max) in self._limits]
        # Scale this axis to the reference
        self._scales = [1.0/self._ranges[i] for i in range(self.n)]
        self._offset = [-min*sca for (min,max), sca in zip(self._limits, self._scales)]

    def _default_ticks(self):
        """Default ticks for added values"""
        ticks = [[min(d),max(d)] for d in numpy.array(self._data).T]
        self.set_ticks(ticks)

    def set_limits(self, limits):
        """Store limits"""
        assert len(limits) == self.n
        self._limits = limits

    def set_ticks(self, values, labels=None):
        """Store ticks"""
        if labels is None: labels=values
        assert len(values) == self.n
        assert len(labels) == self.n
        self._ticks = (values, labels)

    def set_inversion(self, bool):
        """Select which axes whould be inverted"""
        assert len(bool) == self.n
        self._invert = bool

    def _apply_limits(self):
        """Scale the axes to {0,1}"""
        for ax in [self.ax]+self.axes:
            ax.set_ylim(0.0, self._padding_inner + 1.0 + self._padding_outer)

    def _filter_values_and_labels(self, values, labels, filter_large=True, filter_small=True):
        """Filter values and labels - apply floors and ceilings to allowed values"""
        filt_vals = []; filt_labs = []
        filt_upper = self._padding_inner + 1.0
        filt_lower = self._padding_inner + 0.0
        for vals, labs in zip(values, labels):
            vs=[]; ls=[];
            for v, l in zip(vals, labs):
                if (v>filt_lower or not filter_small) and (v<filt_upper or not filter_large):
                    vs.append(v); ls.append(l)
                elif v>=filt_upper:
                    vs.append(filt_upper); ls.append(l)
                elif v<=filt_lower:
                    vs.append(filt_lower); ls.append(l)
            filt_vals.append(vs); filt_labs.append(ls)
        return (filt_vals, filt_labs)

    def _ha_va(self, i_axis, invert=False):

        if (invert is False):
            has = ['center', 'right', 'left']
            vas = ['center', 'bottom', 'top']
        else:
            has = ['center', 'left', 'right']
            vas = ['center', 'top', 'bottom']

        i_frac = float(i_axis) / float(self.n)

        ha = (
            has[0] if (
                ( i_frac <= 0.1 ) or ( i_frac >= 0.9 ) or ( 0.4 <= i_frac <= 0.6 )
            )
            else
            has[1] if (
                ( 0.1 < i_frac < 0.4 )
            )
            else
            has[2]
        )

        va = (
            vas[0] if (
                ( 0.2 <= i_frac <= 0.3 ) or ( 0.7 <= i_frac <= 0.8 )
            )
            else
            vas[1] if (
                ( i_frac < 0.2 ) or ( i_frac > 0.8 )
            )
            else
            vas[2]
        )
        return ha, va

    def _apply_ticks(self):
        tick_vals = self._normalise(self._ticks[0])
        tick_labs = self._ticks[1]
        tick_vals, tick_labs = self._filter_values_and_labels(tick_vals, tick_labs, filter_large=True, filter_small=True)
        for i, (ax, angle, vals, labs) in enumerate(
                zip(
                    self.axes,
                    self.angles,
                    tick_vals,
                    tick_labs,
                )
            ):
            ha, va = self._ha_va(i)
            ax.set_rgrids(
                [v+0.08 for v in vals],
                labels = labs,
                angle = angle,
                fontsize = self.label_fontsize,
                ha = ha,
                va = va,
                bbox = None,
                weight = "bold",
                color = "black",
            )
            # Hide grid lines
            ax.spines["polar"].set_visible(False)

            if (self.fancy_labels is True):
                bbox = dict(boxstyle="round,pad=0.2",facecolor='lightcoral')
                for label in ax.yaxis.get_ticklabels():
                    label.set_bbox(bbox)

    def add(self, values, *args, **kwgs):
        assert len(values) == self.n
        self._data.append(values)
        self._args.append(args)
        self._kwgs.append(kwgs)
        self.d = len(self._data)

    def legend(self):

        self.ax.legend(
            loc = 'upper right',
            bbox_to_anchor = (1.0, 0.0),
            bbox_transform = self.fig.transFigure,
            fontsize = self.legend_fontsize,
            ncol = 2,
        )

    def axis_limits(self):

        assert len(self._limits) == self.n
        assert list(map(len,self._limits)).count(2) == self.n

        lim_vals = self._normalise(self._limits)
        lim_labs = [[round(v,2) for v in vs] for vs in self._limits]

        for i, (vals, labs, ang) in enumerate(
                zip(
                    lim_vals,
                    lim_labs,
                    numpy.deg2rad(numpy.r_[self.angles, self.angles[0]]),
                ),
            ):
            invt = bool(self._invert[i])
            ha, va = 'center', 'center'
            dlta = 0.115 * [1, -1][int(invt)]
            self.ax.text(
                x = ang,
                y = vals[0] - dlta,
                s = labs[0],
                horizontalalignment = ha,
                verticalalignment = va,
                fontsize = self.label_fontsize,
                weight = 'bold',
            )
            self.ax.text(
                x = ang,
                y = vals[1] + dlta,
                s = labs[1],
                horizontalalignment = ha,
                verticalalignment = va,
                fontsize = self.label_fontsize,
                weight = 'bold',
            )

    def show(self):
        if not self._plotted:
            self.plot()
        self.fig.show()

    def plot(self):

        # Calculate defaults where appropriate
        if not self._limits: self._default_limits()
        if not self._scales: self._calculate_scales()
        if not self._ticks:  self._default_ticks()

        # Apply limits, ticks and plot
        self._apply_limits()
        self._apply_ticks()
        self._plotted = True

        norm_vals = numpy.array(
            self._normalise(
                numpy.array(self._data).T.tolist()
            )
        ).T.tolist()

        filt_vals, _ = self._filter_values_and_labels(
            norm_vals,
            norm_vals,
            filter_large = False,
            filter_small = True,
        )

        for values, args, kwgs in zip(filt_vals, self._args, self._kwgs):
            angle = numpy.deg2rad(numpy.r_[self.angles, self.angles[0]])
            vals = numpy.r_[values, values[0]]
            self.ax.plot(angle, vals, *args, **kwgs)

    def savefig(self, filename):
        if not self._plotted:
            self.plot()
        self.fig.savefig(
            filename,
            bbox_inches = "tight",
            pad_inches = 0.5,
        )

    def close(self):
        pyplot.close(self.fig)
