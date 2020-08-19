from matplotlib import pyplot
pyplot.style.use('ggplot')

import numpy

    # List of filenames to return
    filenames = []
    # Get the maximum x-axis
    max_x = max(map(len, plot_vals))
    # Iterate through different images
    for i_img in range(int(1.0*(num_sites-1)/graphs_per_image) + 1):
        # Create filename from template
        # Create filename from template
        f_name = f_template.format(i_img+1)
        filenames.append(f_name)
        # Get new image
        fig, all_axes = pyplot.subplots(graphs_per_image, sharex=False)
        # Change to list if only one subplot (pyplot...)
        if graphs_per_image == 1: all_axes = [all_axes]
        # Iterate through different graphs
        for i_plot, i_site in enumerate(range(graphs_per_image*i_img, graphs_per_image*i_img+graphs_per_image)):
            # Get the axis for this plot
            s_axis = all_axes[i_plot]
            # Change X-axis
            s_axis.set_xticks([-100,1]+range(5,max_x,5)+[max_x])
            s_axis.set_xlim([0,max_x+2])
            # Past the end of the data?
            if i_site >= num_sites:
                s_axis.set_yticks([])
                continue
            # Bar heights
            bar_vals = plot_vals[i_site]
            num_vals = len(bar_vals)
            # Left sides of the bars
            bar_left = [x+0.6 for x in range(num_vals)]
            bar_hght = bar_vals
            # Colour of the bars
            if colour_bool:   bar_colr = ['limegreen' if b else 'red' for b in colour_bool[i_site]]
            elif colour_vals: bar_colr = colour_vals[i_site]
            else:             bar_colr = ['blue']*num_vals
            # Make the bar
            bar_bar = s_axis.bar(left=bar_left, height=bar_hght, width=0.8, color=bar_colr)
            if   colour_bool: s_axis.set_title('Site {}: {} models'.format(i_site+1, sum(colour_bool[i_site])))
            elif colour_vals: s_axis.set_title('Site {}: {} models'.format(i_site+1, colour_vals[i_site].count('limegreen')))
            else:             s_axis.set_title('Site {}'.format(i_site+1))
            # Change Y-axis
            s_axis.set_yticks([int(max(bar_hght)+0.5)])
        # Save
        pyplot.tight_layout()
        pyplot.savefig(f_name)
        pyplot.close(fig)
    return filenames
