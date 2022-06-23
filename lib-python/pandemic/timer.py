import time

class Timer(object):
    def __init__(self):
        self._start = None
        self._last = None

    def get_start_time(self):
        return self._start

    def start(self):
        assert self._start is None, 'Timer already started'
        self._start = time.time()
        return self

    def report_str(self):
        assert self._start is not None, 'Timer must be started first!'
        cur_time = time.time()
        # Report time since start
        out_str = 'Time elapsed since start: {}'.format(self.format_time_delta(cur_time-self._start))
        # Report time since last check
        if self._last is not None:
            out_str += '\nTime elapsed since last check: {}'.format(self.format_time_delta(cur_time-self._last))
        # Store current time as last checkpoint
        self._last = cur_time
        return out_str

    def format_time(self, seconds):
        return time.ctime(seconds)

    def format_time_delta(self, seconds):
        t = int(seconds)
        hours = (t // 3600)
        mins = (t % 3600) // 60
        secs = (t % 60)
        return '{:02d}:{:02d}:{:02d}'.format(hours, mins, secs)
