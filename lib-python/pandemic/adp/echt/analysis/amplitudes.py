


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

