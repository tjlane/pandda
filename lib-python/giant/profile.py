from future import standard_library
standard_library.install_aliases()
import giant.logs as lg
logger = lg.getLogger(__name__)

import cProfile, pstats, io

class profile_code(object):
    def __init__(self):
        self.start()

    def start(self):
        self.profiler = cProfile.Profile()
        self.profiler.enable()

    def stop(self, print_stats=True):
        self.profiler.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(self.profiler, stream=s).sort_stats(sortby)
        if print_stats:
            ps.print_stats()
            logger(s.getvalue())
        return ps
