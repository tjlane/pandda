import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy
import itertools, collections
from libtbx import adopt_init_args

try:
    import matplotlib as mpl
    import matplotlib.cm
    mpl.use('Agg')
    #mpl.interactive(False)
    import matplotlib.pyplot as plt
    plt.switch_backend('agg') # yes I know this is done twice!
    plt.style.use('ggplot')
    #plt.interactive(0)
    plt.rc('font', family='monospace')
except Exception as e:
    logger(e)

def make_linkage_graphs(
    cluster_info,
    output_filename,
    ):

    merging_steps = cluster_info.merging_steps
    n_steps = len(merging_steps)

    if n_steps == 0:
        return False

    linkages = list(zip(*merging_steps))[0]
    labels = list(zip(*merging_steps))[1]

    markers = numpy.arange(0, n_steps, max(1, int(float(n_steps) / 10.)))
    markers[-1] = n_steps-1
    
    fig, axis = plt.subplots(nrows=1, ncols=1)
    axis.set_title('Iterative Clustering')
    axis.plot(linkages, 'b-*', markevery=list(markers))
    axis.set_xlim([-1, n_steps])
    axis.set_ylim([0, 1])
    axis.set_xticks(markers)
    axis.set_xticklabels(markers+1)
    axis.set_xlabel('Clustering Step')
    axis.set_ylabel('Similarity between joined clusters')
    fig.tight_layout()
    fig.savefig(output_filename)
    plt.close(fig)

    return


class Dendrogram(object):


    def __init__(self, 
        n_points,
        children_hash,
        merging_steps,
        inverse_distance = True,
        ):
        _merging_linkages = []
        _merging_pairs = []
        for i, (linkage, pair) in enumerate(merging_steps): 
            #if i > 0: 
            #    # Ensure linkage doesn't decrease
            #    linkage = max(linkage, _merging_linkages[-1])
            _merging_linkages.append(linkage)
            _merging_pairs.append(pair)
        adopt_init_args(self, locals())

    def __call__(self, 
        output_filename,
        i_merge=None,
        ):

        self.x_offset = 1.

        if i_merge == None: 
            i_merge = len(self.merging_steps) - 1
        
        # Extract info for this cluster
        i_cluster = self.n_points + i_merge
        n_children = len(self.children_hash[i_cluster])
        join_pair = self._merging_pairs[i_merge]
        
        root_linkage = self._merging_linkages[i_merge]

        # intial line stub
        x_centre = self.x_offset + float(n_children - 1) / 2.

        fig, axis = plt.subplots(nrows=1, ncols=1)

        if self.inverse_distance is False: 
            self.root_y = 1.05 * root_linkage
            self.leaf_y = 0.0
        else:
            self.root_y = 0.0 # 0.95 * root_linkage
            self.leaf_y = 1.05 * max(self._merging_linkages[:(i_merge+1)])
        
        axis.plot([x_centre, x_centre], [root_linkage, self.root_y])

        # Plot ll line segments below
        for i_join in join_pair: 
            self.iterate_points(
                axis = axis,
                i_cluster = i_join,
                previous_xy = (x_centre, root_linkage),
                )
        
        fig.savefig(output_filename)
        plt.close(fig)

    def iterate_points(self, axis, i_cluster, previous_xy=None):

        prev_x, prev_y = previous_xy

        n_children = len(self.children_hash[i_cluster])

        if n_children == 1: 
            this_x = self.x_offset
            this_y = self.leaf_y
            self.x_offset += 1.
        else: 
            i_merge = i_cluster - self.n_points 
            
            linkage = self._merging_linkages[i_merge]
            join_pair = self._merging_pairs[i_merge]
            
            this_x = self.x_offset + (float(n_children - 1.) /2.)
            this_y = linkage

            #assert not this_y > prev_y
            for i_join in join_pair: 
                self.iterate_points(
                    axis = axis,
                    i_cluster = i_join,
                    previous_xy = (this_x, this_y),
                    )

        # Plot horizontal line
        axis.plot([this_x, prev_x], [prev_y, prev_y], 'b')
        # Plot Vertical line
        axis.plot([this_x, this_x], [this_y, prev_y], 'b')

        return


class LevelPlotAccumulator(object):


    def __init__(self,
            n_groups,
            colour_map_name = 'rainbow',
            x_label = 'Existing Groups',
            y_label = 'Threshold',
            y_width = 1.0,
            y_lim = None,
            title = None,
            fig_axis_pair = None,
            n_colours = 20,
            max_x_labels = 30,
            ):
        adopt_init_args(self, locals())

        self.colour_map = mpl.cm.get_cmap(colour_map_name)

        colours = list(self.get_systematic_colours(n_colours))
        spliced_colours = []
        n_splices = max(1, int(float(n_colours)/5.))
        [spliced_colours.extend(colours[i::n_splices]) for i in range(n_splices)]
        self.colour_cycle = itertools.cycle(spliced_colours)
        self.used_colours = []

        if self.fig_axis_pair is None:
            self.fig, self.axis = plt.subplots(nrows=1, ncols=1)
        else:
            self.fig, self.axis = self.fig_axis_pair

        self.block_counter = 1

        if title is not None: 
            self.axis.set_title(label=str(title))
        
        if (y_lim is not None): 
            y_range = abs(y_lim[1] - y_lim[0]) + y_width
        else: 
            y_range = None

        if n_groups > 100:
            self.font_size = 0
        elif (y_range is not None): 
            self.font_size = 10.0 * (20.0 / y_range)
        else: 
            self.font_size = 10.0

        self.axis.set_xlabel(x_label)
        self.axis.set_ylabel(y_label)

        xticks = list(range(1, n_groups+1))
        self.max_x_labels = min(n_groups, self.max_x_labels)
        plot_every = max(1, 1+((n_groups-1)//self.max_x_labels))
        self.axis.set_xticks(xticks)
        self.axis.set_xticklabels([v if not i%plot_every else '' for i,v in enumerate(xticks)])

        self.axis.set_yticks([])
        self.axis.set_yticklabels([])

        self.legend = None
        self.extra_artists = []

    def get_random_colours(self, n):
        i = numpy.random.random(n)
        return self.colour_map(i)

    def get_systematic_colours(self, n):
        i = numpy.linspace(0., 1., n)
        return self.colour_map(i)        

    def add_line(self, y_value, linestyle='-', linewidth=1, color='k'):

        self.axis.plot((0, self.n_groups+1), (y_value,y_value), linestyle=linestyle, linewidth=linewidth, color=color)

    def add_level(self, list_of_lists_of_indices_tuples, y_value, label=None):

        n_sub_levels = len(list_of_lists_of_indices_tuples)

        bar_width = self.y_width / float(n_sub_levels)

        # Enclosing bar
        self.axis.broken_barh(
            [(0.25, self.n_groups+0.5)], 
            (y_value - self.y_width/2., self.y_width), 
            edgecolor='k', facecolor='w',
            )

        for i_sub_level, level_indices in enumerate(list_of_lists_of_indices_tuples):

            y_sub_level = y_value - (self.y_width/2.) + (i_sub_level * bar_width)

            for i_group, group_indices in enumerate(level_indices):

                n_elements = len(group_indices)

                group_colour = next(self.colour_cycle)
                self.used_colours.append(group_colour)
                
                group_label = self.block_counter

                # Create list of bars (x_min, x_width)
                x_values = [i+1.0 for i in group_indices]
                x_ranges = [(x-0.5, 1.0) for x in x_values]
        
                self.axis.broken_barh(
                    x_ranges, 
                    (y_sub_level, bar_width), 
                    edgecolor='k', facecolor=group_colour,
                    )

                if (self.font_size > 0): #and (n_sub_levels < 3):
                    for i, x in enumerate(x_values): 
                        self.axis.text(
                            x, y_sub_level + bar_width / 2.0, int(group_label), 
                            fontsize=int(float(self.font_size) / float(n_sub_levels)),
                            horizontalalignment='center',
                            verticalalignment='center',
                            rotation=0, rotation_mode="default",
                            )
                
                if False and (self.legend is None): 
                    from matplotlib import patches
                    hdl = patches.Patch(edgecolor='k', facecolor='b')
                    self.legend = self.axis.legend(
                        handles=[hdl], 
                        labels=['New groupings'], 
                        #bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
                        bbox_to_anchor=(0.5, 0.0),
                        bbox_transform=self.fig.transFigure,
                        loc=9, borderaxespad=0.,
                        )
                    self.extra_artists.append(self.legend)

                # Increment
                self.block_counter += 1

        if label is not None: 
            y_ticks = list(self.axis.get_yticks()) + [y_value]
            y_tick_labels = list(self.axis.get_yticklabels()) + [label]
            self.axis.set_yticks(y_ticks)
            self.axis.set_yticklabels(y_tick_labels)

    def simplify_yaxis_labels(self, max_labels=15):

        y_ticks = list(self.axis.get_yticks())
        y_labels = list(self.axis.get_yticklabels())

        if (len(y_ticks) == 0) or (len(y_labels) == 0):
            return

        min_val = min(y_ticks)
        max_val = max(y_ticks)
        max_labels = max(2, max_labels)

        if (len(y_ticks) <= max_labels) or (len(y_labels) <= max_labels):
            return

        min_spacing = (max_val - min_val) / float(max_labels-1)

        y_ticks_out = [y_ticks.pop(0)]
        y_labels_out = [y_labels.pop(0)] 

        last_y = y_ticks_out[0]

        for tick, label in zip(y_ticks, y_labels): 

            if abs(tick - last_y) < min_spacing: 
                continue

            # New group - finish old group
            y_ticks_out.append(tick)
            y_labels_out.append(label)

            last_y = tick
        
        self.axis.set_yticks(y_ticks_out)
        self.axis.set_yticklabels(y_labels_out)

    def format(self):
        self.axis.set_xlim((0, self.n_groups+1))

        for l in self.axis.get_yticklabels():
            l.set_multialignment('right')

        plt.setp(self.axis.get_xticklabels(), rotation=90)

    def write(self, filename):
        assert self.fig is not None
        self.format()
        self.fig.tight_layout()
        self.fig.savefig(
            filename, 
            bbox_extra_artists=self.extra_artists,
            bbox_inches='tight',
            dpi=300)
        plt.close(self.fig)

