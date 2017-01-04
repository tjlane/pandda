import os, sys

import scipy.cluster

from bamboo.common import Info
from giant.stats.cluster import find_connected_groups, generate_group_idxs

from scitbx.array_family import flex

class PointCluster(object):
    def __init__(self, id, points, values):
        """Class to hold information about an identified cluster of points, with associated values"""

        # Check that both are lists for the moment
        points = list(points)
        values = list(values)
        assert len(points) == len(values)

        self.parent = None
        self.id = id
        # Store as flex arrays
        self.points = flex.vec3_double(points)
        self.values = flex.double(values)

        stats = self.values.min_max_mean()
        self.min = stats.min
        self.max = stats.max
        self.mean = stats.mean
        self.size = len(self.values)

        self.peak = self.points[values.index(max(values))]
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
    _attributes = ['estimated_pseudo_occupancy', 'estimated_bdc']
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
            for a in self._attributes: assert hasattr(info, a)
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
    _attributes = [ 'centroid', 'num_events', 'approx_size',
                    'nearest_residues', 'near_crystal_contacts']
    def __init__(self, events=None, id=None, info=None):
        """Class to hold information about an identified site (collection of events)"""
        # List of Events
        self.children = []
        self.parent = None
        self.id = id
        # Add Meta to the object
        if info:
            assert isinstance(info, Info)
            for a in self._attributes: assert hasattr(info, a)
            self.info = info
        else:
            self.info = Info(self._attributes)

        # Add Events
        if events: self.add_events(events=events, update=True)

    def add_events(self, events, update=True):
        if isinstance(events, list):
            assert isinstance(events[0], Event), 'Added events must be of class Event'
            self.children.extend(events)
        else:
            assert isinstance(events, Event), 'Added events must be of class Event'
            self.children.append(events)
        if update: self.update()
        return self

    def apply_parentage(self):
        """Set the site as the parents of the events"""
        for e in self.children: e.parent = self
        return self

    def update(self):
        centroids = flex.vec3_double([e.cluster.centroid for e in self.children])
        self.info.centroid = centroids.mean()
        self.info.approx_size = (centroids.min(), centroids.max())
        self.info.num_events = len(self.children)
        return self

    def summary(self):
        out = []
        out.append('Site Summary')
        out.append('{!s} Events'.format(self.info.num_events))
        out.append('Centroid: {!s}'.format(self.info.centroid))
        out.append('Approx Range: {!s}'.format(self.info.approx_size))
        for s in self.children:
            out.append('-----------------------')
            out.append(s.summary())
        return '\n'.join(out)

    def sort(self, key, reverse=True):
        self.children.sort(key=key, reverse=reverse)
        return self

    def find_protein_context(self, hierarchy):
        return

class SiteList(object):
    def __init__(self, sites=None):
        """Class to hold information about multiple sites on a protein"""
        # List of Sites
        self.children = []
        self.parent = None
        self.id = None
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

def cluster_events(events, cutoff=10, linkage='average'):

    if len(events) == 1:
        return SiteList([Site(events=events, id=1).apply_parentage()])

    centroids = [e.cluster.centroid for e in events]
    cluster_ids = scipy.cluster.hierarchy.fclusterdata( X = centroids,
                                                        t = cutoff,
                                                        criterion = 'distance',
                                                        metric    = 'euclidean',
                                                        method    = linkage )
    cluster_ids = list(cluster_ids)

    sites = []
    for s_idx, e_idxs in generate_group_idxs(cluster_ids):
        assert s_idx > 0
        new_site = Site([events[i] for i in e_idxs], id=s_idx).apply_parentage()
        sites.append(new_site)

    return SiteList(sites)

def event_distance_centroid(event1, event2):
    return (flex.double(event1.cluster.centroid) - flex.double(event2.cluster.centroid)).norm()

def event_distance_min_cluster_distance(event1, event2):
    return min([min((event1.cluster.points-x).norms()) for x in event2.cluster.points])

