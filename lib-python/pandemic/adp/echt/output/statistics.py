import pandas as pd
import pathlib as pl

from scitbx.array_family import flex


class CalculateBFactorStatistics(object):

    def __init__(self):
        pass

    def __call__(self,
        level_b_values,
        level_names,
        structure_factory,
        overall_atom_mask,
        ):

        total_table_dicts = []
        chain_table_dicts = []

        level_b_values = [flex.double(b) for b in level_b_values]

        assert len(level_names) == len(level_b_values)

        total_b_mean = sum([flex.mean(b) for b in level_b_values])
        chain_b_means = {}

        for l_name, l_bvals in zip(level_names, level_b_values):

            hierarchy = structure_factory.custom_copy(
                iso = l_bvals,
                mask = overall_atom_mask,
            ).select(
                overall_atom_mask,
            )

            atom_cache = hierarchy.atom_selection_cache()

            level_b_mean = flex.mean(hierarchy.atoms().extract_b())

            if (total_b_mean > 0.0):
                level_b_mean_perc = 100.0 * level_b_mean / total_b_mean
            else:
                level_b_mean_perc = 0.0

            total_table_dicts.append({
                'Level' : l_name,
                'Chain' : 'all',
                'Average B' : level_b_mean,
                'Average B (%)' : level_b_mean_perc,
            })

            chain_ids = sorted(set([c.id for c in hierarchy.chains()]))

            for c_id in chain_ids:

                chain_sel = atom_cache.selection('chain {}'.format(c_id))
                chain_atoms = hierarchy.atoms().select(chain_sel)
                chain_level_b_mean = flex.mean(chain_atoms.extract_b())

                chain_b_means[c_id] = chain_level_b_mean + chain_b_means.get(c_id, 0)

                chain_table_dicts.append({
                    'Level' : l_name,
                    'Chain' : c_id,
                    'Average B' : chain_level_b_mean,
                })

        for d in chain_table_dicts:

            chain_id = d['Chain']
            chain_b_mean = chain_b_means[chain_id]

            if chain_b_mean > 0.0:
                chain_level_b_mean_perc = 100.0 * d['Average B'] / chain_b_mean
            else:
                chain_level_b_mean_perc = 0.0

            d['Average B (%)'] = chain_level_b_mean_perc

        # Sort in place
        chain_table_dicts.sort(key=lambda d: d['Chain'])

        b_factor_statistics = pd.DataFrame(
            data = (total_table_dicts + chain_table_dicts),
            columns = ['Level', 'Chain', 'Average B', 'Average B (%)'],
        )

        return b_factor_statistics


