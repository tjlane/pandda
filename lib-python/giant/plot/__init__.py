import giant.logs as lg
logger = lg.getLogger(__name__)

from giant.plot.radar import Radar
from giant.plot.simple import *

from giant.exceptions import Sorry

def setup_plotting(no_fail=False, interactive=False, backend='agg', style='ggplot'):

    try:

        import matplotlib
        matplotlib.interactive(interactive)

        from matplotlib import pyplot
        if (pyplot.isinteractive() != interactive):
            logger.warning(
                'Interactive setting is incorrect (current {} is not requested {})'.format(
                    pyplot.isinteractive(),
                    interactive,
                )
            )

        if (backend is not None):
            pyplot.switch_backend(backend)
            if (pyplot.get_backend() != backend):
                logger.warning(
                    'pyplot backend loaded ({}) is not the one requested ({})'.format(
                        pyplot.get_backend(),
                        backend,
                    )
                )

        pyplot.style.use(style)

    except Exception as e:
        if (no_fail is True):
            logger.warning('Failed to setup matplotlib -- may not be able to generate graphs')
        else:
            raise Sorry('Failed to setup matplotlib -- needed for plotting. \nError: {}'.format(str(e)))

    logger.info(
        'matplotlib & pyplot loaded successfully. Using backend "{!s}"'.format(
            pyplot.get_backend(),
        )
    )

