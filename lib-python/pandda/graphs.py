import giant.logs as lg
logger = lg.getLogger(__name__)

try:
    import matplotlib as mpl
    mpl.use('agg')
    mpl.interactive(False)
    from matplotlib import pyplot as plt
    plt.switch_backend('agg') # yes I know this done twice -- for safety!
    plt.interactive(0)
except:
    logger(str(e))
    import matplotlib as mpl
    from matplotlib import pyplot as plt

import json, collections
import pathlib as pl
import numpy as np

def failure_graph(title):
    fig = plt.figure()
    plt.title('Failed to make {}'.format(title))
    return fig

class PlotConfig(object):

    def __init__(self):

        self.colour_map_name = None
        self.plot_style = None
        self.font_family = None
        self.font_name = None

    def configure(self,
        colour_map_name = None,
        plot_style = None,
        font_family = None,
        font_name = None,
        ):

        if (colour_map_name is not None):
            self.colour_map_name = colour_map_name

        if (plot_style is not None):
            self.set_plot_style(
                plot_style = plot_style,
                )

        if (font_family is not None):
            self.set_font_family(
                font_family = font_family,
                )

        if (font_name is not None):
            self.set_font(
                font_name = font_name,
                font_family = font_family,
                )

    def set_plot_style(self,
        plot_style,
        ):

        if (plot_style == 'xkcd'):

            try:

                plt.xkcd()

                logger('Theme set to "xkcd".')

            except Exception as e:

                logger.warning('Failed to set plot style "xkcd".')

            self.plot_style = plot_style

            return None

        try:

            logger(
                'Setting plot_style to "{}"'.format(plot_style)
                )

            plt.style.use(plot_style)

            self.plot_style = plot_style

        except Exception as e:

            logger.warning('Failed to set plot style to {}.\n\t{}'.format(plot_style, str(e)))

    def set_font_family(self,
        font_family,
        ):

        try:

            logger(
                'Setting font_family to "{}"'.format(font_family)
                )

            plt.rc('font', family=font_family)
            
            self.font_family = font_family

        except Exception as e:

            logger.warning(
                'Failed to set font family to {}.\n\t{}'.format(font_family, str(e))
                )

    def set_font(self,
        font_name,
        font_family,
        ):

        try:

            if font_family is None:
                logger.warning(
                    'Cannot set font: must provide a font family in order to set font.'
                    )
                return None

            font_family_str = 'font.'+str(font_family)

            if font_family_str not in list(plt.rcParams.keys()):
                logger.warning(
                    'Cannot set font: invalid font family provided "{}".'.format(font_family)
                    )
                return None

            family_fonts = plt.rcParams[font_family_str]

            if (font_name not in family_fonts):
                logger.warning(
                    'font "{}" does not exist in font family "{}". Setting the font may not work. (valid options: {})'.format(
                        font_name, 
                        font_family, 
                        ', '.join(family_fonts),
                        )
                    )

            logger('Setting font_name to "{}"'.format(font_name))

            plt.rcParams[font_family_str].insert(0, font_name)

            self.font_name = font_name

        except Exception as e:

            logger.warning('Failed to set font to {}.\n\t{}'.format(font_name, str(e)))

    def get_colour_map(self):
        return mpl.cm.get_cmap(self.colour_map_name)

    def get_colours(self, n):
        cm = self.get_colour_map()
        return cm(np.linspace(0., 1., n))


config = PlotConfig()

def configure(
    colour_map_name = 'rainbow',
    plot_style = 'ggplot',
    font_family = 'monospace',
    font_name = None,
    ):

    global config

    config.configure(
        colour_map_name = colour_map_name,
        plot_style = plot_style,
        font_family = font_family,
        font_name = font_name,
        )

    return config


#####


class JsonPlotter(object):
    """Format graph data into json for html plotting"""
    
    def __init__(self):
        self.data = {}

    def __str__(self):
        return self.as_json()

    def as_json(self):
        return json.dumps(self.data)

    def as_javascript(self):
        return ""


class PanddaPlotter(object):

    output_key = None

    def __init__(self, output_path):

        self.output_path = output_path

    def __call__(self, *args, **kwargs):

        fig = self.plot(*args, **kwargs)

        # try:
        #     fig = self.plot(*args, **kwargs)
        # except:
        #     fig = failure_graph(title='error')

        filename = self.get_path()

        self.save(
            fig = fig,
            filename = filename, 
            )

        return {self.output_key : filename}

    def setup(self, nrows=1, ncols=1, **kwargs):

        fig, axes = plt.subplots(
            nrows=nrows, ncols=ncols,
            **kwargs
            )

        return fig, axes

    def save(self, fig, filename):

        fig.tight_layout()
        fig.savefig(filename)
        plt.close(fig) 

    def get_path(self):

        p = pl.Path(self.output_path)

        if not p.parent.exists():
            p.parent.mkdir(parents=True)

        return str(p)

    def json(self):
        return self.json_plotter.as_json()


class PanddaDatasetPlotter(PanddaPlotter):

    output_key = None

    def __init__(self, output_path_template):
        
        self.output_path_template = str(output_path_template)
        self.json_plotter = JsonPlotter()

        assert ('{label}' in output_path_template) # remove this 

    def __call__(self, datasets, *args, **kwargs):
        
        output_files = collections.OrderedDict()

        for dkey, dataset in sorted(datasets.items()):

            fig = self.plot(
                dataset = dataset,
                dataset_label = dkey,
                *args, **kwargs
                )
             
            # try: 
            #     fig = self.plot(
            #         dataset = dataset,
            #         dataset_label = dkey,
            #         *args, **kwargs
            #         )
            # except:
            #     fig = failure_graph(title='error')

            filename = self.get_path(
                dataset_label = dkey,
                )
            
            self.save(
                fig = fig, 
                filename = filename,
                )

            output_files[dkey] = filename

        return {self.output_key : output_files}

    def get_label(self, dataset=None):

        if (dataset is None): 
            return None

        return dataset.tag

    def get_path(self, dataset_label=None, **kwargs):

        p = pl.Path(
            self.output_path_template.format(
                label = dataset_label,
                **kwargs
                )
            )
        
        if not p.parent.exists():
            p.parent.mkdir(parents=True)

        return str(p)


class PanddaMultiDatasetPlotter(PanddaPlotter):

    output_key = None

    def __init__(self, output_path):
        
        self.output_path = str(output_path)
        self.json_plotter = JsonPlotter()

    def __call__(self, datasets, *args, **kwargs):

        dataset_keys = sorted(datasets.keys())
        dataset_list = [datasets[k] for k in dataset_keys]

        fig = self.plot(
            dataset_labels = dataset_keys,
            dataset_list = dataset_list,
            *args, **kwargs
            )

        # try:
        #     fig = self.plot(
        #         dataset_labels = dataset_keys,
        #         dataset_list = dataset_list,
        #         *args, **kwargs
        #         )
        # except:
        #     fig = failure_graph(title='error')

        filename = self.get_path()

        self.save(
            fig = fig,
            filename = filename,
            )

        return {self.output_key : filename}
