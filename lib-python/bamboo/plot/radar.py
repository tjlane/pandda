from matplotlib import pyplot
pyplot.style.use('ggplot')

import numpy

class Radar(object):

    def __init__(self, titles, rect=None, fig=None):
        # Make a rectangle
        if rect is None: rect = [0.05, 0.05, 0.95, 0.95]
        # Number of plots
        self.n = len(titles)
        # Number of data points
        self.d = 0

        # Control flags
        self._plotted    = False

        # Inside and Outside offsets
        self._fixed_offset = 0.2
        self._fixed_buffer = 0.2

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
        if fig: self.fig = fig
        else:   self.fig = pyplot.figure()

        # Store titles
        self.titles = titles
        # Create and store angles and axes
        self.angles = [a if a <=360. else a - 360. for a in numpy.arange(90, 90+360, 360.0/self.n)]
        self.axes   = [self.fig.add_axes(rect, projection="polar", label="axes%d" % i) for i in range(self.n)]

        # Get main axis, add theta grids
        self.ax = self.axes[0]
        self.ax.set_thetagrids(self.angles, labels=titles, frac=1.2, fontsize=14, weight="bold", color="black")
        # Remove all things from other axes
        for ax in self.axes[1:]:
            ax.patch.set_visible(False)
            ax.grid("off")
            ax.xaxis.set_visible(False)
            self.ax.yaxis.grid(False)

    def _normalise(self, vals):
        """Normalise values to the scaled axis"""
        assert len(vals) == self.n
        ret_vals = []
        for vls,s,o,i in zip(vals, self._scales, self._offset, self._invert):
            ret_vals_sub = []
            for v in vls:
                if (v is None) or numpy.isnan(v): ret_vals_sub.append(self._fixed_offset+0.5)
                else:                             ret_vals_sub.append(self._fixed_offset+(1.0-i)*(v*s+o)+i*(1.0-(v*s+o)))
            ret_vals.append(ret_vals_sub)
        return ret_vals

    def _default_limits(self):
        """Set the limits to the range of the data"""
        limits = []
        for d in numpy.array(self._data).T:
            if   len(d)==1:      limits.append((d[0]-1,d[0]+1))
            elif min(d)==max(d): limits.append((d[0]-1,d[0]+1))
            else:                limits.append((min(d),max(d)))
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
        norm_limits = self._normalise(self._limits)
        for ax, (min, max) in zip(self.axes, norm_limits):
            ax.set_ylim(0.0, 1.0+self._fixed_offset+self._fixed_buffer)

    def _filter_values_and_labels(self, values, labels, filter_large=True, filter_small=True):
        """Filter values and labels - apply floors and ceilings to allowed values"""
        filt_vals = []; filt_labs = []
        filt_upper = self._fixed_offset + 1.0
        filt_lower = self._fixed_offset + 0.0
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

    def _apply_ticks(self):
        tick_vals = self._normalise(self._ticks[0])
        tick_labs = self._ticks[1]
        tick_vals, tick_labs = self._filter_values_and_labels(tick_vals, tick_labs, filter_large=True, filter_small=True)
        for ax, angle, vals, labs in zip(self.axes, self.angles, tick_vals, tick_labs):
            ax.set_rgrids(vals, labels=labs, angle=angle, fontsize=14, weight="bold", color="black")
            ax.spines["polar"].set_visible(False)
            ax.xaxis.grid(True,color='black',linestyle='-')

    def add(self, values, *args, **kwgs):
        assert len(values) == self.n
        self._data.append(values)
        self._args.append(args)
        self._kwgs.append(kwgs)
        self.d = len(self._data)

    def show(self):
        if not self._plotted: self.plot()
        self.fig.show()

    def plot(self):
        if not self._limits: self._default_limits()
        if not self._scales: self._calculate_scales()
        if not self._ticks:  self._default_ticks()
        self._apply_limits()
        self._apply_ticks()
        self._plotted = True
        norm_vals = numpy.array(self._normalise(numpy.array(self._data).T.tolist())).T.tolist()
        filt_vals, dummy = self._filter_values_and_labels(norm_vals, norm_vals, filter_large=False, filter_small=True)
        for values, args, kwgs in zip(filt_vals, self._args, self._kwgs):
            angle = numpy.deg2rad(numpy.r_[self.angles, self.angles[0]])
            vals = numpy.r_[values, values[0]]
            self.ax.plot(angle, vals, *args, **kwgs)

    def savefig(self, filename):
        if not self._plotted: self.plot()
        self.fig.savefig(filename, bbox_inches="tight", pad_inches=0.5)

    def close(self):
        pyplot.close(self.fig)
