import os, sys

from Bamboo.Common import Info

class PointCluster(object):
    def __init__(self, id, points, values):
        """Class to hold information about an identified cluster of points, with associated values"""
        assert len(points) == len(values)
        self.id = id
        self.points = flex.vec3_double(points)
        self.values = flex.double(values)
        self.stats = self.values.min_max_mean()
        self.size = len(self.values)
        self.centroid = self.points.mean()

class Event(object):
    def __init__(self, id, cluster, info=None):
        """Class to hold information about an event in a dataset"""
        self.id = id
        self.cluster = cluster
        if info:
            assert isinstance(info, Info)
            assert hasattr(info, 'estimated_occupancy')
            self.info = info
        else:
            self.info = Info(['estimated_occupancy'])

class Site(object):
    def __init__(self, events=[]):
        """Class to hold information about an identified site (collection of events)"""
        self.children = []
        if events: self.add_events(events)
        self.update()
    def add_events(self, events, update=True):
        if isinstance(events, list):
            assert isinstance(events[0], Event), 'Added events must be of class Event'
            self.children.extend(events)
        else:
            assert isinstance(events, Event), 'Added events must be of class Event'
            self.children.append(events)
        if update:
            self.update()
        return self
    def update(self)
        self.size = len(self.children)
        self.centroid = flex.vec3_double([e.cluster.centroid for e in self.children]).mean()
        return self

class SiteList(object):
    def __init__(self, sites=[]):
        """Class to hold information about multiple sites on a protein"""
        assert isinstance(sites, list)
        if sites:
            assert isinstance(sites[0], Site)
        self.sites = sites
        self.update()
    def add_site(self, site)
        assert isinstance(site, Site), 'Added sites must be of class Site'
        self.sites.append(site)
        self.update()
        return self
    def update(self):
        return self

def cluster_events(events, site_list=None):
    """Cluster events into sites. Returns a SiteList object. Appends new sites to existing SiteList if given"""

    if not site_list: site_list = SiteList()

    for new_event in events:
        create_new_site = True
        if site_list.sites:
            # Check to see if this event should be put in this site
            event_centroid =

        if create_new_site:
            # Create new site and add to site_list
            new_site = Site(new_event)
            site_list.add_site(new_site)


