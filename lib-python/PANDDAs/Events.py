import os, sys

import libtbx.cluster

from Bamboo.Common import Info
from Giant.Stats.Cluster import find_connected_groups

from scitbx.array_family import flex

class PointCluster(object):
    def __init__(self, id, points, values):
        """Class to hold information about an identified cluster of points, with associated values"""
        assert len(points) == len(values)
        self.parent = None
        self.id = id
        self.points = flex.vec3_double(points)
        self.values = flex.double(values)

        stats = self.values.min_max_mean()
        self.min = stats.min
        self.max = stats.max
        self.mean = stats.mean

        self.size = len(self.values)

        self.peak = self.points[self.values.index(max(self.values))]
        self.centroid = self.points.mean()

    def summary(self):
        out = []
        out.append('Point Cluster {!s} Summary'.format(self.id))
        out.append('{!s} Points'.format(self.size))
        out.append('Min:'.format(self.min))
        out.append('Max:'.format(self.max))
        out.append('Mean:'.format(self.mean))
        return '\n'.join(out)

class Event(object):
    _attributes = ['estimated_occupancy']
    def __init__(self, id, cluster, info=None):
        """Class to hold information about an event in a dataset"""
        # Process cluster object (points and values for Event)
        assert isinstance(cluster, PointCluster), 'cluster must be of type PointCluster'
        self.cluster = cluster
        self.cluster.parent = self
        # Give it a name
        if id: self.id = id
        else:  self.id = cluster.id
        # Allow a parent
        self.parent = None
        # Add Meta to the object
        if info:
            assert isinstance(info, Info)
            for a in self._attributes:
                assert hasattr(info, a)
            self.info = info
        else:
            self.info = Info(self._attributes)

    def summary(self):
        out = []
        out.append('Event {!s} Summary'.format(self.id))
        out.append(self.cluster.summary())
        out.append(self.info.summary())
        return '\n'.join(out)

class Site(object):
    def __init__(self, events=None):
        """Class to hold information about an identified site (collection of events)"""
        # List of Events
        self.children = []
        self.parent = None
        self.id = None
        # Add Events
        if events: self.add_events(events)

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

    def update(self):
        # Number of events at the site
        self.num_events = len(self.children)
        # Extract site centroids
        centroids = flex.vec3_double([e.cluster.centroid for e in self.children])
        # Centroid of the site (mean of event centroids)
        self.centre = centroids.mean()
        # Size of the site (min->max)
        self.approx_size = (centroids.min(), centroids.max())
        return self

    def summary(self):
        out = []
        out.append('Site Summary')
        out.append('{!s} Events'.format(self.num_events))
        out.append('Centroid: {!s}'.format(self.centre))
        out.append('Approx Range: {!s}'.format(self.approx_size))
        for s in self.children:
            out.append('-----------------------')
            out.append(s.summary())
        return '\n'.join(out)

    def sort(self, key, reverse=True):
        return self.children.sort(key=key, reverse=True)

class SiteList(object):
    def __init__(self, sites=None):
        """Class to hold information about multiple sites on a protein"""
        # List of Sites
        self.children = []
        # Add Sites
        if sites: self.add_sites(sites)

    def add_sites(self, sites, update=True):
        if isinstance(sites, list):
            assert isinstance(sites[0], Site), 'Added sites must be of class Site'
            self.children.extend(sites)
        else:
            assert isinstance(sites, Site), 'Added sites must be of class Site'
            self.children.append(sites)
        if update:
            self.update()
        return self

    def update(self):
        # Number of events at the site
        self.num_sites = len(self.children)
        return self

    def summary(self):
        out = []
        out.append('SiteList Summary')
        out.append('{!s} Sites'.format(self.num_sites))
        for s in self.children:
            out.append('=====================================')
            out.append(s.summary())
        return '\n'.join(out)

    def sort(self, key, reverse=True):
        self.children.sort(key=key, reverse=True)
        return self

    def renumber(self, start_at=1):
        for i_c, c in enumerate(self.children): c.id=i_c+start_at
        return self

def event_distance(event1, event2):
    return (flex.double(event1.cluster.centroid) - flex.double(event2.cluster.centroid)).norm()

def cluster_events(events, cutoff=5, linkage='average'):
    """Cluster events into sites. Returns a SiteList object. Appends new sites to existing SiteList if given"""

    # Create hierarchical clustering of the events
    cl = libtbx.cluster.HierarchicalClustering(data=events, distance_function=event_distance, linkage=linkage)
    # Group the events
    clusters = cl.getlevel(cutoff)
    # Turn lists of events into sites
    sites = map(Site, clusters)
    for i_site, site in enumerate(sites):
        site.id = i_site+1
        for event in site.children:
            event.parent = site
    # Merge into site list
    site_list = SiteList(sites)
    return site_list




