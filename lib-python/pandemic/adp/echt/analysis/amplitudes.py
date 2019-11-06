amplitude_directory = easy_directory(os.path.join(output_directory, 'group_amplitude_distributions'))


class AnalyseAmplitudes:

    

    def calculate_amplitudes_dendrogram(self, out_dir_tag):
        """Cluster the amplitudes across the datasets"""

        fm = self.file_manager
        out_dir = fm.get_dir(out_dir_tag)
        dendro_template = fm.add_file(file_name='amplitude_dendrogram-level_{}.png', file_tag='png-amplitude-dendrogram-template', dir_tag=out_dir_tag)

        self.log.heading('Generating dendrograms of TLS amplitudes between datasets')

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
            self.log('For level {}:'.format(level.index))
            self.log('\tNumber of groups: {}'.format(n_grp))
            self.log('\tNumber of TLS: {}'.format(n_tls))
            self.log('\tNumber of datasets: {}'.format(n_dst))

            amp_mat = group_amplitudes.reshape((n_grp*n_tls, n_dst)).T
            all_amplitudes.append(amp_mat)

            # Plot dendrogram
            fname = dendro_template.format(level.index)
            self.log('> {}'.format(fname))
            link_mat = scipy.cluster.hierarchy.linkage(amp_mat, method='average', metric='euclidean')
            dendrogram(fname=fname, link_mat=link_mat, labels=labels)

        # Plot dendrogram
        fname=dendro_template.format('all')
        self.log('> {}'.format(fname))
        link_mat = scipy.cluster.hierarchy.linkage(numpy.concatenate(all_amplitudes, axis=1), method='average', metric='euclidean')
        dendrogram(fname=fname, link_mat=link_mat, labels=labels)

    def show_group_amplitude_statistics(self,
        group_amplitudes,
        ):

        log = self.log

        log.subheading('Summary of group amplitudes for levels in ECHT model')
        for i_level, ga in enumerate(group_amplitudes):
            log.bar()
            log('Amplitude statistics by dataset: Level {}'.format(i_level+1))
            log.bar(False, True)
            # Plot histogram of amplitudes for each dataset for this level
            ga_new = ga.copy().set_index(['group','mode','label'])
            ga_new.index = range(len(ga_new.index))
            ga_new = ga_new.transpose()
            for label, values in ga_new.iterrows():
                log.bar()
                # Plot histogram of amplitudes across all groups
                if len(values) > 5:
                    try:
                        from ascii_graph import Pyasciigraph
                        g=Pyasciigraph(float_format='{0:d}')
                        counts, bounds = numpy.histogram(a=values, bins=10)
                        graph_data = [('{:.3f}-{:.3f}'.format(bounds[i],bounds[i+1]), v) for i,v in enumerate(counts)]
                        for l in g.graph(label='\n> Histogram of amplitudes for dataset {}\n'.format(label), data=graph_data):
                            if l.startswith('#######'): continue
                            log(l.replace(u"\u2588", '=').replace('= ','> '))
                    except ImportError:
                        pass
                    except:
                        raise
                # Write summary statistics
                log('\n> Amplitude statistics for level {}, dataset {}\n'.format(i_level+1, label))
                mn = values.mean()
                md = values.median()
                sd = values.std()
                log('> Mean:   {:.3f} A^2 (B-factor {:.3f} A^2)'.format(mn, EIGHT_PI_SQ*mn))
                log('> Median: {:.3f} A^2 (B-factor {:.3f} A^2)'.format(md, EIGHT_PI_SQ*md))
                log('> Std:    {:.3f} A^2 (B-factor {:.3f} A^2)'.format(sd, EIGHT_PI_SQ*sd))
            log.bar(True, True)