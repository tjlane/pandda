import collections
import numpy as np

from giant.common.clustering import (
    generate_group_idxs,
    )


class AssignEventsToSites(object):

    def __init__(self,
        clustering_distance = 5.0,
        clustering_linkage = "average",
        ):

        self.clustering_distance = clustering_distance
        self.clustering_linkage = clustering_linkage

    def __call__(self,
        dataset_event_dicts,
        ):

        if len(dataset_event_dicts) == 0:
            return []

        sites, site_idxs = self.group_events_into_sites(
            events = dataset_event_dicts,
            )

        for i, e in enumerate(dataset_event_dicts):

            # Get the site_idx for this event
            site_idx = site_idxs[i]

            # Get the site
            s = sites[site_idx]

            assert s['site_idx'] == site_idx

            e['site_idx'] = s['site_idx']
            e['site_num'] = s['site_num']

        return sites

    def group_events_into_sites(self, events):

        import scipy.cluster

        if len(events) == 1:
            s = self.make_site(events=events)
            s['site_idx'] = 0
            s['site_num'] = 1
            return [s], [0]

        centroids = [
            e['xyz_centroid_ref']
            for e in events
            ]

        site_ids = list(
            scipy.cluster.hierarchy.fclusterdata(
                X = centroids,
                t = self.clustering_distance,
                criterion = 'distance',
                metric = 'euclidean',
                method = self.clustering_linkage,
                )
            )

        sites = []
        event_site_idxs = np.zeros(len(events), dtype=int)

        for site_idx, (s_id, e_idxs) in enumerate(generate_group_idxs(site_ids)):

            new_site = self.make_site(
                events = [events[i] for i in e_idxs],
                )
            sites.append(new_site)

            event_site_idxs[e_idxs] = site_idx

        # Renumber the sites
        site_list_sort_hash = sorted(
            range(len(sites)),
            key = lambda i: (
                sites[i]['n_events'],
                sites[i]['max_value'],
                ),
            reverse = True,
            )

        # List of sites
        sorted_sites = [
            sites[i]
            for i in site_list_sort_hash
            ]

        # Map events to sites
        event_site_idxs = [
            site_list_sort_hash.index(i)
            for i in event_site_idxs
            ]

        # Assign the sites labels
        for i, s in enumerate(sorted_sites):
            s['site_idx'] = i
            s['site_num'] = i+1

        return sorted_sites, event_site_idxs

    def make_site(self, events):

        site = collections.OrderedDict([
            (
                'site_idx', None,
                ),
            (
                'site_num', None,
                ),
            (
                'n_events', len(events),
                ),
            (
                'max_value', np.max([e['z_peak'] for e in events]),
                ),
            (
                'xyz_centroid', tuple(
                    round(v, 3) for v in
                    np.mean([e['xyz_centroid_ref'] for e in events], axis=0)
                    ),
                ),
            (
                'xyz_extent', (
                    tuple(
                        round(v, 3) for v in
                        np.min([e['xyz_centroid_ref'] for e in events], axis=0)
                        ),
                    tuple(
                        round(v, 3) for v in
                        np.max([e['xyz_centroid_ref'] for e in events], axis=0)
                        ),
                    ),
                ),
            ])

        return site
