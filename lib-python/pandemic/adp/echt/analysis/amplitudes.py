import giant.logs as lg
logger = lg.getLogger(__name__)

from libtbx import adopt_init_args, group_args
from giant.paths import easy_directory


class AnalyseTLSAmplitudesTask(object):


    def __init__(self,
        output_directory,
        ):

        adopt_init_args(self, locals())

    def run(self):
        pass

    def calculate_amplitudes_dendrogram(self, out_dir_tag):
        """Cluster the amplitudes across the datasets"""

        fm = self.file_manager
        out_dir = fm.get_dir(out_dir_tag)
        dendro_template = fm.add_file(file_name='amplitude_dendrogram-level_{}.png', file_tag='png-amplitude-dendrogram-template', dir_tag=out_dir_tag)

        logger.heading('Generating dendrograms of TLS amplitudes between datasets')

        # List of linkage-distances for each level
        level_distance_matrices = []

        labels = [m.tag for m in self.models]

        all_amplitudes = []
        for level in self.fitter.levels:
            # List of amplitudes for each group in this level
            group_amplitudes = []
            for i_group, sel, fitter in level:
                tls_mats, tls_amps = fitter.result()
                group_amplitudes.append(tls_amps)
            group_amplitudes = numpy.array(group_amplitudes)

            n_grp, n_tls, n_dst = group_amplitudes.shape
            logger('For level {}:'.format(level.index))
            logger('\tNumber of groups: {}'.format(n_grp))
            logger('\tNumber of TLS: {}'.format(n_tls))
            logger('\tNumber of datasets: {}'.format(n_dst))

            amp_mat = group_amplitudes.reshape((n_grp*n_tls, n_dst)).T
            all_amplitudes.append(amp_mat)

            # Plot dendrogram
            fname = dendro_template.format(level.index)
            logger('> {}'.format(fname))
            link_mat = scipy.cluster.hierarchy.linkage(amp_mat, method='average', metric='euclidean')
            dendrogram(fname=fname, link_mat=link_mat, labels=labels)

        # Plot dendrogram
        fname=dendro_template.format('all')
        logger('> {}'.format(fname))
        link_mat = scipy.cluster.hierarchy.linkage(numpy.concatenate(all_amplitudes, axis=1), method='average', metric='euclidean')
        dendrogram(fname=fname, link_mat=link_mat, labels=labels)

    def show_group_amplitude_statistics(self,
        group_amplitudes,
        ):

        logger.subheading('Summary of group amplitudes for levels in ECHT model')
        for i_level, ga in enumerate(group_amplitudes):
            logger.bar()
            logger('Amplitude statistics by dataset: Level {}'.format(i_level+1))
            logger.bar(False, True)
            # Plot histogram of amplitudes for each dataset for this level
            ga_new = ga.copy().set_index(['group','mode','label'])
            ga_new.index = list(range(len(ga_new.index)))
            ga_new = ga_new.transpose()
            for label, values in ga_new.iterrows():
                logger.bar()
                # Plot histogram of amplitudes across all groups
                if len(values) > 5:
                    try:
                        from ascii_graph import Pyasciigraph
                        g=Pyasciigraph(float_format='{0:d}')
                        counts, bounds = numpy.histogram(a=values, bins=10)
                        graph_data = [('{:.3f}-{:.3f}'.format(bounds[i],bounds[i+1]), v) for i,v in enumerate(counts)]
                        for l in g.graph(label='\n> Histogram of amplitudes for dataset {}\n'.format(label), data=graph_data):
                            if l.startswith('#######'): continue
                            logger(l.replace(u"\u2588", '=').replace('= ','> '))
                    except ImportError:
                        pass
                    except:
                        raise
                # Write summary statistics
                logger('\n> Amplitude statistics for level {}, dataset {}\n'.format(i_level+1, label))
                mn = values.mean()
                md = values.median()
                sd = values.std()
                logger('> Mean:   {:.3f} A^2 (B-factor {:.3f} A^2)'.format(mn, constants.EIGHTPISQ*mn))
                logger('> Median: {:.3f} A^2 (B-factor {:.3f} A^2)'.format(md, constants.EIGHTPISQ*md))
                logger('> Std:    {:.3f} A^2 (B-factor {:.3f} A^2)'.format(sd, constants.EIGHTPISQ*sd))
            logger.bar(True, True)

                # Write histograms of amplitudes -- only for non-zero models
                #if (tls_matrices.sum() > 0.0) and self.distribution_images:
                #    hist_png = self.filepath(self.level_tls_amplitudes_png.format(i_level+1, i_group+1), self.amplitude_directory)
                #    titles = ['Mode {}:'.format(i_m+1) for i_m in xrange(model_object.n_modes)]
                #    x_vals = [tls_amplitudes[i_m,:]    for i_m in xrange(model_object.n_modes)]
                #    self.plot.multi_histogram(filename  = hist_png,
                #                              x_vals    = x_vals,
                #                              titles    = titles,
                #                              x_labs    = ['']*model_object.n_modes,
                #                              rotate_x_labels = True,
                #                              shape     = (tls_amplitudes.shape[0], 1),
                #                              n_bins    = 30, x_lim=[0, None])
                #    file_dict.setdefault('level_tls_amplitudes_png',collections.OrderedDict())[(level_name,group_object.label)] = amp_csv
