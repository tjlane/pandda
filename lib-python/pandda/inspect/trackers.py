import numpy as np


class iTracker(object):

    def __init__(self,
        n_total,
        i_start = 0,
        ):

        self.i_current = i_start
        self.n_total = n_total

    def get(self):

        return self.i_current

    def set(self, i):

        self.i_current = (
            i % self.n_total
            )

        return self.get()

    def next(self):

        self.set(self.i_current + 1)

        return self.get()

    def prev(self):

        self.set(self.i_current - 1)

        return self.get()

    def at_first(self):

        return self.get() == 0

    def at_last(self):

        return self.get() == (self.n_total - 1)


class EventTracker(iTracker):
    pass


class LigandTracker(iTracker):
    pass


class SiteTracker(object):

    def __init__(self,
        site_idxs,
        event_tracker,
        ):

        self.site_idxs = site_idxs
        self.n_total = len(set(site_idxs))
        
        if not np.isnan(self.site_idxs).all():
            assert self.site_idxs.min() == 0
            assert self.site_idxs.max() == (self.n_total - 1)

        self.event_tracker = event_tracker

    def get(self):

        i_event = self.event_tracker.get()

        return self.site_idxs[i_event]

    def next(self):

        i_event = self.event_tracker.get()

        curr_site_idx = self.site_idxs[i_event]

        next_site_idx = (curr_site_idx + 1) % self.n_total

        i_event_new = self.find_first_event(next_site_idx)

        self.event_tracker.set(i_event_new)

        return next_site_idx

    def prev(self):

        i_event = self.event_tracker.get()

        curr_site_idx = self.site_idxs[i_event]

        prev_site_idx = (curr_site_idx - 1) % self.n_total

        i_event_new = self.find_first_event(prev_site_idx)

        self.event_tracker.set(i_event_new)

        return prev_site_idx

    def find_first_event(self, site_idx):

        event_idxs = np.where(self.site_idxs == site_idx)[0]

        return event_idxs[0]


