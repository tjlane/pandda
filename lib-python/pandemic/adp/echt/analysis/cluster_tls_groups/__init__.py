import giant.logs as lg
logger = lg.getLogger(__name__)

import collections

from libtbx import adopt_init_args, group_args
from libtbx.utils import Sorry

from pandemic.adp.utils import show_file_dict
from pandemic.adp.echt.analysis.cluster_tls_groups.helliger import ClusterTLSGroups_HelligerDistance
from pandemic.adp.echt.analysis.cluster_tls_groups.overlap import ClusterTLSGroups_OverlapMass


class ClusterTLSGroupsTask(object):


    show_file_dict = show_file_dict

    max_groups = 50
    characters_to_print = 90

    _cluster_classes = {
        'helliger_distance' : ClusterTLSGroups_HelligerDistance,
        'overlap_mass'  : ClusterTLSGroups_OverlapMass,
        }

    def __init__(self,
        output_directory,
        parameters,
        metric = 'overlap_mass',
        write_levels_function = None,
        ):
        adopt_init_args(self, locals(), exclude=('parameters',))

        # Extract method-specific class
        try:
            cluster_class = self._cluster_classes[metric]
        except KeyError:
            raise Sorry('Method does not exist: {}'.format(metric))

        # Extract method-specific parameters -- a little hacky at the moment
        if not isinstance(parameters, dict):
            parameters = vars(parameters)
        try:
            metric_parameters = parameters[metric]
        except:
            raise Failure('Failed to extract metric parameters ({}) from input parameters'.format(metric))

        if not isinstance(metric_parameters, dict):
            metric_parameters = vars(metric_parameters)
            metric_parameters = {k:v for k,v in metric_parameters.items() if not k.startswith('_')}

        self.cluster = cluster_class(
            xyz_cutoff_distance = parameters['xyz_cutoff_distance'],
            **metric_parameters
            )

    def as_html_summary(self):
        from pandemic.adp.echt.analysis.cluster_tls_groups.html import ClusterTLSGroupsTaskHtmlSummary
        return ClusterTLSGroupsTaskHtmlSummary(self)

    def run(self,
        model_object,
        model_files,
        ):

        logger.subheading('Clustering TLS Groups to generate new levels & groupings', spacer=True)

        if not self.output_directory.exists():
            self.output_directory.mkdir(parents=True)

        # Extract structures here so that no knowledge of the dictionary key is required at higher levels
        model_structures = model_files['level_uijs_pdb']

        # Output files and objects
        clustering_results = collections.OrderedDict()
        output_files = collections.OrderedDict()

        # Iterate through levels
        for level_name in model_object.tls_level_names:

            # Want absolute index of level
            i_l = model_object.all_level_names.index(level_name)
            i_tls = model_object.tls_level_names.index(level_name)

            # Extract tls group objects from model
            tls_groups = model_object.tls_objects[i_tls]
            n_groups = len(tls_groups)

            # Skip if processing more groups than maximum
            if (n_groups == 1) or (n_groups > self.max_groups):
                continue

            logger.subheading('Clustering groups for level {} ({})'.format(i_l+1, level_name))

            # Perform clustering
            c_result = self.cluster(
                tls_groups = tls_groups,
                )

            # Write clustering output
            c_files = self.cluster.write_output(
                result = c_result,
                tls_groups = tls_groups,
                model_structure = model_structures[level_name],
                output_prefix = 'tls_clustering_level_{}'.format(i_l+1),
                output_directory = self.output_directory,
                )

            # Store output files/objects
            clustering_results[level_name] = c_result
            output_files[level_name] = c_files

        if self.write_levels_function is not None:
            of = self.write_new_levels(
                clustering_results = clustering_results,
                model_object = model_object,
                write_levels_function = self.write_levels_function,
                )
            for k, v in of.items():
                output_files.setdefault(k, {})['new_levels'] = v

        self.show_file_dict(output_files)

        self.result = group_args(
            clustering_results = clustering_results,
            output_files = output_files,
            )

        return self.result

    def write_new_levels(self,
        clustering_results,
        model_object,
        write_levels_function,
        ):

        output_files = collections.OrderedDict()

        out_directory = (self.output_directory / 'new_groupings')

        if not out_directory.exists():
            out_directory.mkdir(parents=True)

        logger.subheading('Outputting levels generated at each clustering threshold', spacer=True)

        for level_name, level_data in clustering_results.items():

            o_dict = output_files.setdefault(level_name, collections.OrderedDict())

            # Want absolute index of level
            i_l = model_object.all_level_names.index(level_name)

            logger.subheading('New groupings derived from level {} ({})'.format(i_l+1, level_name))

            for threshold, hierarchy in level_data.module_info.threshold_unique_hierarchies.items():

                filename = str(out_directory / 'level_{}_threshold_{}.eff'.format(i_l+1, threshold))

                label_template = '{} (groups {{}})'.format(level_name)

                phil_str = write_levels_function(
                    level_name = level_name,
                    indices_hierarchy = hierarchy,
                    output_filename = filename,
                    label_template = label_template,
                    comment_lines = [],
                    depth = None,
                    insert_before = level_name,
                    insert_after = None,
                    reverse_levels = False,
                    )
                o_dict[threshold] = filename

                # Report
                logger.bar()
                logger('@ threshold: {}'.format(threshold))
                logger.bar()
                logger(self.truncate_lines(phil_str, max_lines=10))
                #logger(phil_str[:self.characters_to_print] + '...'*(len(phil_str)>self.characters_to_print))
                logger('> Written to {}'.format(filename))

        return output_files

    def truncate_lines(self, string, max_lines=None):
        lines = string.split('\n')
        output_lines = [l[:self.characters_to_print]+'...'*(len(l)>self.characters_to_print) for l in lines]
        if (max_lines is not None) and (len(output_lines) > max_lines):
            output_lines = output_lines[:max_lines] + ['...']
        return '\n'.join(output_lines)





