import giant.logs as lg
logger = lg.getLogger(__name__)

import collections 
import numpy as np
import pandas as pd


class WriteEchtModelTables(object):

    level_tls_matrices_csv   = 'tls_matrices_level_{level_num:04d}.csv'
    level_tls_amplitudes_csv = 'tls_amplitudes_level_{level_num:04d}.csv'
    level_tls_origins_csv    = 'tls_origins_level_{level_num:04d}.csv'

    def __init__(self, 
        output_directory,
        ):

        self.output_directory = output_directory


    def __call__(self,
        model_object,
        ):
    
        if not self.output_directory.exists():
            self.output_directory.mkdir(parents=True)

        output_files = collections.OrderedDict()

        output_files['level_tls_matrices_csv'] = self.write_matrices_table(
            model_object = model_object,
            )
        
        output_files['level_tls_amplitudes_csv'] = self.write_amplitudes_table(
            model_object = model_object,
            )
        
        output_files['level_tls_origins_csv'] = self.write_origins_table(
            model_object = model_object,
            )

        return output_files

    def write_matrices_table(self, model_object):

        output_files = collections.OrderedDict()

        for i_level, tls_objects in enumerate(model_object.tls_objects):
            
            level_name = model_object.tls_level_names[i_level]

            table_rows = []

            for i_group, group_object in enumerate(tls_objects):

                tls_matrices = np.array(
                    [mode.matrices.get() for mode in group_object.tls_parameters]
                    )

                assert tls_matrices.shape == (model_object.n_modes, 21)

                for i_mode in range(model_object.n_modes):

                    headers = np.array(
                        [ i_group+1, i_mode+1, group_object.label ],
                        dtype = object,
                        )
                    
                    table_rows.append(
                        np.concatenate(
                            [ headers, tls_matrices[i_mode] ]
                            )
                        )

            table = pd.DataFrame(
                columns = [
                    "group", "mode", "label",
                    "T11","T22","T33","T12","T13","T23",
                    "L11","L22","L33","L12","L13","L23",
                    "S11","S12","S13","S21","S22","S23","S31","S32","S33",
                    ],
                data = table_rows,
                # dtype = object,
                )

            path = str(
                self.output_directory / self.level_tls_matrices_csv.format(
                    level_num = i_level+1,
                    )
                )

            table.to_csv(path)

            output_files[level_name] = path

        return output_files

    def write_amplitudes_table(self, model_object):
        
        output_files = collections.OrderedDict()

        for i_level, tls_objects in enumerate(model_object.tls_objects):

            level_name = model_object.tls_level_names[i_level]

            table_rows = []

            for i_group, group_object in enumerate(tls_objects):

                tls_amplitudes = np.array(
                    [mode.amplitudes.get() for mode in group_object.tls_parameters]
                    )

                assert tls_amplitudes.shape == (model_object.n_modes, model_object.n_datasets)

                for i_mode in range(model_object.n_modes):

                    headers = np.array(
                        [ i_group+1, i_mode+1, group_object.label ], 
                        dtype = object,
                        )

                    table_rows.append(
                        np.concatenate(
                            [ headers, tls_amplitudes[i_mode] ]
                            )
                        )

            table = pd.DataFrame(
                columns = ["group", "mode", "label"] + list(model_object.dataset_labels), 
                data = table_rows,
                dtype = object,
                )

            path = str(
                self.output_directory / self.level_tls_amplitudes_csv.format(
                    level_num = i_level+1,
                    )
                )

            table.to_csv(path)

            output_files[level_name] = path

        return output_files
    
    def write_origins_table(self, model_object):

        output_files = collections.OrderedDict()
        
        for i_level, tls_objects in enumerate(model_object.tls_objects):

            level_name = model_object.tls_level_names[i_level]

            table_rows = []

            for i_group, group_object in enumerate(tls_objects):

                tls_origins = np.array(
                    list(map(str, group_object.origins.round(3))),
                    dtype = str,
                    )

                assert tls_origins.shape == (model_object.n_datasets,)

                headers = np.array(
                    [ i_group+1, group_object.label ], 
                    dtype = object,
                    )

                table_rows.append(
                    np.concatenate(
                        [ headers, tls_origins ]
                        )
                    )

            table = pd.DataFrame(
                columns = ["group", "label"] + list(model_object.dataset_labels), 
                data = table_rows,
                dtype = object,
                )

            path = str(
                self.output_directory / self.level_tls_origins_csv.format(
                    level_num = i_level+1,
                    )
                )

            table.to_csv(path)

            output_files[level_name] = path

        return output_files

    
    
