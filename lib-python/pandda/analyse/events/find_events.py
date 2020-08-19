import collections 


class BasicPanddaFindEvents: 

    sort_key = "cluster_size"

    def __init__(self, 
        z_map_type,
        find_clusters,
        filter_clusters,
        ):

        self.z_map_type = z_map_type
        self.find_clusters = find_clusters
        self.filter_clusters = filter_clusters

    def __call__(self,
        map_grid,
        map_mask,
        search_mask,
        dataset_key,
        dataset,
        dataset_map_data,
        statistical_model,
        ):

        z_map_data = statistical_model.get_z_map(
            dataset_key = dataset_key,
            map_data = dataset_map_data,
            z_map_type = self.z_map_type,
            )

        # Index the analysis mask on the map mask
        reindex_mask, reindex_sel = search_mask.index_on_other(
            other = map_mask,
            )

        # Embed
        z_map = reindex_mask.embed_data(
            data = z_map_data[reindex_sel],
            )

        clusters = self.find_clusters(
            map_data = z_map,
            )

        clusters = self.filter_clusters(
            clusters = clusters,
            map_grid = map_grid,
            dataset = dataset,
            )

        events = [
            self.make_event(
                cluster = cluster,
                map_grid = map_grid,
                dataset_key = dataset_key,
                dataset = dataset,
                ) 
            for cluster in clusters
            ]

        # Sort events
        events = sorted(
            events,
            key = lambda e: e[self.sort_key],
            )

        # Number the output events
        for i, e in enumerate(events):
            e['event_idx'] = i
            e['event_num'] = i+1

        return events

    def make_event(self, 
        cluster,
        map_grid,
        dataset_key,
        dataset,
        ):

        event = collections.OrderedDict([
            (
                'dtag', dataset_key,
                ),
            (
                'event_idx', None, # populated later
                ),
            (
                'event_num', None, # populated later
                ),
            (
                'site_idx', None, # populated later
                ),
            (
                'site_num', None, # populated later
                ),
            (
                'event_fraction', None, # populated later
                ), 
            (
                'bdc', None, # populated later
                ), 
            (
                'z_peak', round(cluster.max(),3),
                ),
            (
                'z_mean', round(cluster.mean(),3),
                ),
            (
                'cluster', cluster,
                ),
            (
                'cluster_size', cluster.size(),
                ),
            (
                'xyz_peak', tuple(
                    round(v, 3) for v in
                    dataset.model.alignment.ref2nat(
                        map_grid.grid2cart([cluster.peak()])
                        )[0]
                    ),
                ),
            (
                'xyz_centroid', tuple(
                    round(v, 3) for v in
                    dataset.model.alignment.ref2nat(
                        map_grid.grid2cart([cluster.centroid()])
                        )[0]
                    ),
                ),
            (
                'xyz_peak_ref', tuple(
                    round(v, 3) for v in
                    map_grid.grid2cart([cluster.peak()])[0]
                    ),
                ),
            (
                'xyz_centroid_ref', tuple(
                    round(v, 3) for v in
                    map_grid.grid2cart([cluster.centroid()])[0]
                    ),
                ),
            ])
                
        return event
