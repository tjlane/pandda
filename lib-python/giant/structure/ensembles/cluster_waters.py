import giant.logs as lg
logger = lg.getLogger(__name__)

from giant.common import EIGHTPISQ

import giant.common.geometry as gm

import os
import numpy as np
from scitbx.array_family import flex

from giant.processors import ProcessorJoblib
import multiprocessing


from matplotlib import pyplot as plt


class GetWaters(object):

    def __init__(self):

        pass

    def __call__(self, hierarchy):

        ac = hierarchy.atom_selection_cache()

        h = hierarchy.select(
            ac.selection('resname HOH')
            )

        return h


class FilterWaters(object):

    def __init__(self):

        pass

    def __call__(self, hierarchy):

        # from scitbx.array_family import flex

        # sel = (hierarchy.atoms().extract_b() < 0.)

        # i_sel = flex.size_t(
        #     np.random.randint(
        #         0,
        #         sel.size(),
        #         10000,
        #         )
        #     )

        # sel.set_selected(i_sel, True)

        # hierarchy = hierarchy.select(sel)

        return hierarchy


class AtomCluster(object):

    def __init__(self, points):

        self.points = np.array(points)

    def __len__(self):

        return self.size()

    def __iter__(self):

        for p in self.points:
            yield p

    def distance_to(self, point):

        return gm.distance_between(point, self.centroid())

    def minimum_distance_to(self, points):

        pair_dists = gm.pairwise_distances(self.points, points)

        return pair_dists.min()

    def centroid(self):

        return self.points.mean(axis=0)

    def extent(self):

        p_dists = gm.pairwise_distances(self.points)

        max_dist = p_dists.max()

        return max_dist

    def msf(self):

        rmsf = self.rmsf()

        return (rmsf * rmsf)

    def rmsf(self):

        deltas = (
            self.points - self.centroid()
            )

        rmsf = gm.rms_coordinates(deltas)

        return rmsf

    def covariance(self):

        deltas = (
            self.points - self.centroid()
            )

        x,y,z = deltas.T

        xx = np.mean(x*x)
        yy = np.mean(y*y)
        zz = np.mean(z*z)

        xy = np.mean(x*y)
        yz = np.mean(y*z)
        zx = np.mean(z*x)

        covariance = np.array([
            [xx, xy, zx],
            [xy, yy, yz],
            [zx, yz, zz],
            ])

        return covariance

    def b_iso(self):

        msf = self.msf()

        return (
            EIGHTPISQ * msf
            )

    def u_iso(self):

        return self.msf()

    def uij_iso(self):

        msf = self.msf()

        return (
            msf,
            msf,
            msf,
            0.0,
            0.0,
            0.0,
            )

    def uij_aniso(self):

        cov = self.covariance()

        return (
            cov[0,0],
            cov[1,1],
            cov[2,2],
            cov[0,1],
            cov[0,2],
            cov[1,2],
            )

    def size(self):

        return self.points.shape[0]

    def as_atom(self):

        from iotbx.pdb.hierarchy import atom

        a = atom()
        a.set_b(self.b_iso())
        a.set_uij(self.uij_aniso())
        a.set_xyz(self.centroid())

        return a


class AtomClusterList(object):

    def __init__(self, atom_clusters):

        self.atom_clusters = atom_clusters

    def __iter__(self):

        for c in self.atom_clusters:
            yield c

    def centroids(self):

        return np.array([c.centroid() for c in self])

    def extents(self):

        return np.array([c.extent() for c in self])

    def msfs(self):

        return np.array([c.msf() for c in self])

    def rmsfs(self):

        return np.array([c.rmsf() for c in self])

    def covariances(self):

        return np.array([c.covariance() for c in self])

    def b_isos(self):

        return np.array([c.b_iso() for c in self])

    def u_isos(self):

        return np.array([c.u_iso() for c in self])

    def uij_isos(self):

        return np.array([c.uij_iso() for c in self])

    def uij_anisos(self):

        return np.array([c.uij_aniso() for c in self])

    def sizes(self):

        return np.array([c.size() for c in self])

    def as_atoms(self):

        return np.array([c.as_atom() for c in self])


class ClusteringResult(object):

    def __init__(self,
        best_model,
        all_models,
        ):

        self.best_model = best_model
        self.all_models = all_models


class ClusterPoints(object):

    def __init__(self,
        remove_singletons = True,
        ):

        self.remove_singletons = bool(
            remove_singletons
            )

    def __call__(self, points):

        points = np.array(points)

        clustering_result = self.cluster(
            points = points,
            )

        selections = self.assign_points(
            points = points,
            fitted_model = clustering_result.best_model,
            )

        # logger(
        #     'Clustered {} points into {} clusters'.format(
        #         len(points),
        #         len(selections),
        #         )
        #     )

        return (selections, clustering_result)

    def __str__(self):

        return self.name

    def assign_points(self, points, fitted_model):

        labels = self.predict(
            points = points,
            fitted_model = fitted_model,
            )

        set_labels = set(labels)

        if self.remove_singletons is True:
            set_labels.difference_update([-1])

        selections = [
            (labels == i)
            for i in
            sorted(set_labels)
            ]

        return selections


class ClusterPointsDBScan(ClusterPoints):

    name = "ClusterPointsDBScan"

    def __init__(self,
        remove_singletons = True,
        dbscan_eps = 1.0,
        dbscan_min_samples = 6, # 2*dim
        ):

        self.remove_singletons = bool(remove_singletons)
        self.dbscan_eps = float(dbscan_eps)
        self.dbscan_min_samples = float(dbscan_min_samples)

    def cluster(self, points):

        dbscan = self.fit(
            points = points,
            )

        return ClusteringResult(
            best_model = dbscan,
            all_models = None,
            )

    def fit(self, points):

        from sklearn.cluster import DBSCAN

        dbscan = DBSCAN(
            eps = self.dbscan_eps,
            min_samples = self.dbscan_min_samples,
            )

        dbscan.fit(points)

        return dbscan

    def predict(self, points, fitted_model):

        return fitted_model.labels_


class ClusterPointsGMM(ClusterPoints):

    name = "ClusterPointsGMM"

    def __init__(self,
        remove_singletons = True,
        n_atoms_max = 50,
        n_atoms_steps = (10, 5, 2, 1),
        gaussian_type = "spherical",
        ):

        self.remove_singletons = bool(remove_singletons)

        self.n_atoms_max = int(n_atoms_max)
        self.n_atoms_steps = tuple(map(int,n_atoms_steps))

        self.gaussian_type = gaussian_type

        self.validate()

    def validate(self):

        assert self.gaussian_type in ["spherical", "tied", "full"]

    def get_n_range(self, min_value, max_value, step_size, n_min=1, n_max=None):

        if n_min is not None:
            min_value = max(1, min_value)

        if n_max is not None:
            max_value = min(n_max, max_value)

        n_range = [min_value, max_value] + list(np.arange(
            (min_value // step_size) * step_size + step_size,
            max_value,
            step_size,
            ))

        return sorted(set(n_range))

    def cluster(self, points):

        n_points = len(points)

        n_max = min(n_points, self.n_atoms_max)

        n_limits = (
            1,
            n_max,
            )

        all_gmms = []

        for n_step in self.n_atoms_steps:

            if n_step > n_points:
                continue

            n_values = self.get_n_range(
                min_value = n_limits[0],
                max_value = n_limits[1],
                step_size = n_step,
                n_min = 1,
                n_max = n_max,
                )

            # logger(
            #     'Limits: {}'.format(
            #         str(tuple(n_values))
            #         )
            #     )

            best_gmm, all_gmms_step = self.optimised_fit(
                points = points,
                n_values = n_values,
                )

            n_limits = (
                best_gmm.n_components - n_step,
                best_gmm.n_components + n_step,
                )

            all_gmms.extend(all_gmms_step)

        all_gmms = sorted(
            all_gmms,
            key = lambda g: g.n_components,
            )

        return ClusteringResult(
            best_model = best_gmm,
            all_models = all_gmms,
            )

    def fit(self, points, n_components):

        from sklearn import mixture

        gmm = mixture.GaussianMixture(
            n_components = n_components,
            covariance_type = self.gaussian_type,
            )

        gmm.fit(points)

        return gmm

    def predict(self, points, fitted_model):

        return fitted_model.predict(points)

    def optimised_fit(self, points, n_values):

        best_gmm = None
        best_bic = np.infty

        all_gmms = []

        for n_components in n_values:

            if n_components > points.shape[0]:
                continue

            gmm = self.fit(
                points = points,
                n_components = n_components,
                )

            all_gmms.append(gmm)

            bic = gmm.bic(points)

            is_better = False

            if bic < best_bic:
                best_gmm = gmm
                best_bic = bic
                is_better = True

            logger.debug(
                'GMM with {} components:  BIC: {:>20}   | {:>20}   {}'.format(
                    n_components,
                    bic,
                    np.log(bic),
                    '*'*int(is_better),
                    )
                )

        return (best_gmm, all_gmms)


class AtomClusterer(object):

    def __init__(self, cluster_func):

        self.cluster_func = cluster_func

    def __call__(self, points):

        if len(points) == 1:
            # Singleton
            return [ [True] ]

        cluster_selections, cluster_result = self.cluster_func(
            points,
            )

        return cluster_selections


class WriteOutput(object):

    def __init__(self):

        pass

    def __call__(self, input_hierarchy, atom_clusters, out_dir):

        if out_dir is not None:
            if not out_dir.exists():
                out_dir.mkdir(parents=True)

        atom_clusters = sorted(
            atom_clusters,
            key=lambda c: c.size(),
            reverse=True,
            )

        # currently assumes cctbx shared_atom
        atom_cluster_list = AtomClusterList([
            AtomCluster(
                points = a.extract_xyz(),
                )
            for a in atom_clusters
            ])

        self.pymol_script(
            input_hierarchy = input_hierarchy,
            atom_cluster_list = atom_cluster_list,
            out_path = out_dir / 'pymol_script.py',
            )

        self.cluster_histograms(
            input_hierarchy = input_hierarchy,
            atom_cluster_list = atom_cluster_list,
            out_path = out_dir / 'histograms.png',
            )

    def pymol_script(self, input_hierarchy, atom_cluster_list, out_path):

        from giant.pymol_utils import (
            PythonScript, shapes
            )

        n_models = input_hierarchy.models_size()

        p = PythonScript()

        for i, atom_clust in enumerate(atom_cluster_list):

            occ = min(
                1.0,
                float(len(atom_clust)) / float(n_models),
                )

            color = (
                min(1.0, occ),
                0.,
                max(0.0, 1.0-occ),
                )

            obj = 'cluster{}'.format(i+1)

            p.add_shapes(
                cgo_shapes = [
                    shapes.Sphere(
                        centre = pos,
                        radius = 0.4,
                        color = color,
                        )
                    for pos in atom_clust
                    ],
                obj = obj,
                )

            if occ < 0.1:
                p.disable(obj=obj)

        # Make COMs
        for i, atom_clust in enumerate(atom_cluster_list):

            p.add_generic_object(
                item = shapes.Pseudoatom(
                    pos = atom_clust.centroid(),
                    vdw = 0.1,
                    label = str(i+1),
                    color = (1.,1.,1.),
                    ),
                obj = 'cluster{}_com'.format(i+1),
                )

            # p.label(
            #     selection='cluster{}_com'.format(i+1),
            #     expression=str(i+1),
            #     )

        p.set("label_size", 25)
        p.set("label_position", (0,0,4))
        p.set("label_color", "white", "all")

        if out_path.exists():
            os.remove(str(out_path))

        p.write_script(str(out_path))

    def cluster_histograms(self, input_hierarchy, atom_cluster_list, out_path):

        n_models = input_hierarchy.models_size()

        fig, axes = plt.subplots(2,2)

        (a1,a2,a3,a4) = axes.flatten()

        a1.hist(
            atom_cluster_list.sizes(),
            bins = 30,
            )
        a1.set_xlabel('Cluster size')
        a1.set_ylabel('Count')

        a2.scatter(
            atom_cluster_list.sizes(),
            atom_cluster_list.extents(),
            )
        a2.set_xlabel('SIZE')
        a2.set_ylabel('Extents')

        a3.hist(
            atom_cluster_list.rmsfs(),
            bins = 30,
            )
        a3.set_xlabel('RMSF')
        a3.set_ylabel('Count')

        a4.scatter(
            atom_cluster_list.sizes(),
            atom_cluster_list.rmsfs(),
            )
        a4.set_xlabel('SIZE')
        a4.set_ylabel('RMSF')

        a4.set_title(
            'Correlation: {}'.format(
                np.corrcoef([
                        atom_cluster_list.sizes(),
                        atom_cluster_list.rmsfs(),
                    ])[0,1]
                ),
            )

        plt.tight_layout()
        plt.savefig(str(out_path), dpi=200)
        plt.close(fig)


class ClusterWaters(object):

    def __init__(self,
        cluster_steps = None,
        min_occupancy = 0.01,
        ):

        self.get_waters = GetWaters()

        self.filter_waters = FilterWaters()

        self.cluster_steps = (
            cluster_steps
            if
            (cluster_steps is not None)
            else
                [
                ClusterPointsDBScan(),
                ClusterPointsGMM(),
                ]
            )

        self.write_output = WriteOutput()

        self.min_occupancy = float(min_occupancy)

        self.processor = ProcessorJoblib(
            n_cpus = multiprocessing.cpu_count(),
            )

    def __call__(self, hierarchy, out_dir=None):

        n_states = hierarchy.models_size()

        min_atoms_per_cluster = int(
            max(
                1,
                np.floor(self.min_occupancy * float(n_states))
                )
            )

        logger('Extracting waters')

        waters_h = self.get_waters(hierarchy)

        logger('Filtering waters')

        hierarchy = self.filter_waters(waters_h)

        logger('Clustering atoms')

        atom_clusters = [waters_h.atoms()]

        for i_step, cluster_points in enumerate(self.cluster_steps):

            logger.subheading(
                'Step {}: {}'.format(
                    i_step+1,
                    str(cluster_points),
                    )
                )

            output_atom_clusters = self.apply_clustering(
                cluster_func = cluster_points,
                atom_clusters = atom_clusters,
                )

            if out_dir is not None:

                self.write_output(
                    input_hierarchy = hierarchy,
                    atom_clusters = output_atom_clusters,
                    out_dir = (
                        out_dir / 'step{}'.format(i_step+1)
                        ),
                    )

            logger(
                'Input clusters: {}\nOutput clusters: {}'.format(
                    len(atom_clusters),
                    len(output_atom_clusters),
                    )
                )

            atom_clusters = self.filter_clusters(
                atom_clusters = output_atom_clusters,
                min_atoms_per_cluster = min_atoms_per_cluster,
                )

        return atom_clusters

    def apply_clustering(self, cluster_func, atom_clusters):

        clusterer = AtomClusterer(
            cluster_func = cluster_func,
            )

        clustering_selections = self.processor([
            self.processor.make_wrapper(
                func = clusterer,
                points = a.extract_xyz(),
                )
            for a in atom_clusters
            ])

        new_atom_clusters = []

        for atom_selections, atoms in zip(clustering_selections, atom_clusters):

            for a_sel in atom_selections:
                new_atom_clusters.append(
                    atoms.select(flex.bool(a_sel))
                    )

        return new_atom_clusters

    def filter_clusters(self, atom_clusters, min_atoms_per_cluster=1):

        input_n = len(atom_clusters)

        atom_clusters = [
            a
            for a in atom_clusters
            if (len(a) > min_atoms_per_cluster)
            ]

        output_n = len(atom_clusters)

        logger(
            'filtered {} clusters to {} clusters based on minimum size of {} points'.format(
                input_n, output_n, min_atoms_per_cluster,
                )
            )

        return atom_clusters
