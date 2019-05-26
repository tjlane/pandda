from libtbx import adopt_init_args
from bamboo.common.logs import Log


class WarningLogger:


    def __init__(self, log=None):
        _warnings = []
        # Counter to record the last warning printed
        _i = 0 
        if log is None: log = Log()
        adopt_init_args(self, locals())

    def __call__(self, warnings):
        if isinstance(warnings, str): 
            self.append(warnings)
        else:
            self.extend(warnings)
        self.flush()

    def append(self, warning, show=False):
        self._warnings.append(warning)
        if show: self.flush()

    def extend(self, warnings, show=False):
        self._warnings.extend(warnings)
        if show: self.flush()

    def flush(self):
        """Print unshown messages"""
        _i = self._i
        n = len(self._warnings) - _i
        self.log.subheading('{} Warnings'.format(n))
        if n == 0: 
            return
        for w in self._warnings[_i:]:
            self.log.bar(True, True)
            self.log(w)
        self.log.bar(True, True)
        self._i = n

    def report(self):
        """Print all messages"""
        n = len(self._warnings)
        if n == 0:
            self.log.subheading('Completed with 0 warnings')
            return
        banner = '{} non-fatal errors/warnings occurred (see below)'.format(n)
        self.log.subheading(banner)
        for i, e in enumerate(self._warnings):
            self.log.bar()
            self.log('Message {} of {}'.format(i+1, n))
            self.log.bar()
            self.log(e)
        self.log.bar()
        self.log.subheading(banner.replace('below','above'))