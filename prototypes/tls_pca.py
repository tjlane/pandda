class MultiDatasetTLSAnalysis(object):

    _t_name = ['T11','T22','T33','T12','T13','T23']
    _l_name = ['L11','L22','L33','L12','L13','L23']
    _s_name = ['S11','S12','S13','S21','S22','S23','S31','S32','S33']
    _all_names = _t_name+_l_name+_s_name

    pca_groups = [('T',   _t_name),
                  ('L',   _l_name),
                  ('S',   _s_name),
                  ('TLS', _all_names)]

    def __init__(self, models, out_dir='tls-analysis', csv_base='tls-params-'):

        self.out_dir = easy_directory(out_dir)

        self.csv_base = os.path.join(self.out_dir, os.path.basename(csv_base))
        self.tables = {}

        self.log = Log(verbose=True)

        for m in models:
            self.add(m)

    def add(self, model):

        tls_params = extract_tls_from_pdb(model.filename)

        for tls_fit in tls_params.tls_params:
            # Extract data table for this selection
            tls_table = self.tables.setdefault(tls_fit.selection_string, pandas.DataFrame(columns=self._t_name+self._l_name+self._s_name))

            tls_table.loc[model.tag]=None
            tls_table.loc[model.tag,self._t_name] = tls_fit.t
            tls_table.loc[model.tag,self._l_name] = tls_fit.l
            tls_table.loc[model.tag,self._s_name] = tls_fit.s

    def show(self):
        for selection in sorted(self.tables.keys()):
            self.log.subheading(selection)
            print(self.tables[selection])

    def write(self):
        with open(self.csv_base+'selections.log', 'w') as csv_log:
            for i, selection in enumerate(sorted(self.tables.keys())):
                n = i+1
                csv_log.write('selection {:03d} : {}\n'.format(n, selection))
                self.tables[selection].to_csv(self.csv_base+'selection-{:03d}'.format(n)+'.csv')

    def make_plots(self):

        from IPython import embed; embed()

        raise SystemExit()

    def run_pca(self):

        # Iterate through the different TLS selections
        for tls_sel_str in sorted(self.tables.keys()):

            sel_table = self.tables[tls_sel_str]

            # Analyse different parts of the TLS parameterisation
            for group_name, group_sel in self.pca_groups:

                self.log.subheading('{}   -   ({})'.format(group_name, ','.join(group_sel)))

                sel_data = sel_table[group_sel].values.astype('float64')

                pca = mdp.nodes.PCANode(reduce=True)
                pca.train(sel_data)
                pca.stop_training()

                # Bolt on a couple of other statistics
                pca.d_perc = 100.0*pca.d/pca.total_variance
                pca.component_corr_to_avg = numpy.corrcoef(pca.avg, pca.get_recmatrix())[0][1:]

                # Summary
                self.log('The PCA was trained on {} {} matrices (of {} parameters each)'.format(pca.tlen, group_name, pca.input_dim))
                self.log('The resultant PCA is formed of {} principal components'.format(pca.output_dim))
                self.log('...which explain {}% of the observed variance'.format(100.0*pca.explained_variance))
                self.log('...and which individually explain')
                self.log('\t{}'.format('\n\t'.join(['{:7.3f}\t({:7.3f}%)'.format(v,p) for v,p in zip(pca.d,pca.d_perc)])))
                self.log('...of the variance.')
                self.log.bar()
                self.log('The average TLS parameters are')
                self.log('\t'+'\n\t'.join(['{} : {:8.3f}'.format(*v) for v in zip(group_sel,pca.avg.tolist()[0])]))
                self.log.bar()
                self.log('The correlations between the average TLS model and each of the principle components is:')
                self.log('\t{}'.format('\n\t'.join(['{:8.3f}'.format(s) for s in pca.component_corr_to_avg])))

        #        print 'projection matrix: '
        #        print pca.get_projmatrix()

                proj_vals = pca.execute(sel_data)

                print(proj_vals.shape)

        #        self.plot_bar(pca.d)
        #        self.plot_2d(proj_vals[:,:2])

                from IPython import embed; embed()

                #help(pca)

    def plot_2d(self, data, block=True):
        fig = pyplot.figure(figsize=(8,8))
        ax = fig.add_subplot(111)
        pyplot.rcParams['legend.fontsize'] = 10
        ax.plot(data[:,0], data[:,1], 'o', markersize=8, color='blue', alpha=0.5)
        ax.set_aspect('equal')
        pyplot.show(block)

    def plot_3d(self, data, block=True):
        fig = pyplot.figure(figsize=(8,8))
        ax = fig.add_subplot(111, projection='3d')
        pyplot.rcParams['legend.fontsize'] = 10
        ax.plot(data[:,0], data[:,1], data[:,2], 'o', markersize=8, color='blue', alpha=0.5)
        ax.set_aspect('equal')
        pyplot.show(block)

    def plot_bar(self, vals, block=True):
        fig = pyplot.figure(figsize=(8,8))
        ax = fig.add_subplot(111)
        pyplot.rcParams['legend.fontsize'] = 10
        ax.plot(range(1,len(vals)+1), vals, 'o', markersize=8, color='blue', alpha=0.5)
        pyplot.show(block)

