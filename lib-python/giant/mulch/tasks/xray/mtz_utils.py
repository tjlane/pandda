import giant.logs as lg
logger = lg.getLogger(__name__)

import os, shutil, traceback

import pathlib as pl

from giant.mulch.tasks.utils import (
    run_program,
    raise_missing,
    )

from giant.mulch.tasks.io import (
    TaskReturn,
    TaskReturnStatus,
    ModelDataInputOutput,
    )

from giant.mulch.labelling import (
    PathLabeller,
    )

from giant.dispatcher import (
    Dispatcher,
    )

from giant.mulch.dataset import (
    CrystallographicData,
    )


class RemoveColumns(object):

    def __init__(self, 
        columns_to_remove,
        out_suffix = ".remove_columns.mtz",
        ):

        self.columns_to_remove = columns_to_remove
        self.out_suffix = out_suffix

    def __call__(self, 
        in_mtz,
        ):

        logger.subheading('Removing columns')

        out_mtz = in_mtz.with_suffix(self.out_suffix)

        prog = Dispatcher('cad')

        prog.extend_args([
            'hklin1', str(in_mtz),
            'hklout', str(out_mtz),
        ])

        prog.extend_stdin([
            'monitor BRIEF',
        ])

        input_columns = CrystallographicData.from_file(
            str(in_mtz),
            ).crystal.column_labels

        keep_columns = [
            col
            for col in input_columns
            if col not in (
                ['H','K','L'] + self.columns_to_remove
                )
        ]

        keep_comms = [
            "E{i}={col}".format(
                i = i+1,
                col = col
                )
            for i, col 
            in enumerate(keep_columns)
        ]

        prog.extend_stdin([
            'labin file_number 1 {com}'.format(
                com = ' '.join(keep_comms),
                ),
            'labout file_number 1 {com}'.format(
                com = ' '.join(keep_comms),
                ),
        ])

        run_program(prog)

        if not out_mtz.exists():

            prog.write_output(
                str(out_mtz.with_suffix('.log'))
                )
            
            raise_missing(
                filepath = out_mtz, 
                result = prog.result,
                )

        return out_mtz


class ReindexMtzToReference(object):
    """
    Reindex the data in one mtz to a reference mtz
    """

    def __init__(self, 
        reference_mtz,
        out_suffix = ".reindex.mtz",
        tolerance = 5,
        ):

        self.reference_mtz = reference_mtz
        self.out_suffix = out_suffix
        self.tolerance = tolerance

        assert self.reference_mtz is not None

    def __call__(self, 
        in_mtz, 
        ):

        logger.subheading('Running reindexing')

        out_mtz = in_mtz.with_suffix(self.out_suffix)

        prog = Dispatcher('pointless')

        prog.extend_args([
            'hklin',  str(in_mtz),
            'hklref', str(self.reference_mtz),
            'hklout', str(out_mtz),
        ])

        prog.extend_stdin([
            'tolerance {}'.format(self.tolerance),
        ])

        run_program(prog)

        if not out_mtz.exists():

            prog.write_output(
                str(out_mtz.with_suffix('.log'))
                )

            raise_missing(
                filepath = out_mtz,
                result = prog.result,
                )

        return out_mtz


class FillMissingReflections(object): 

    def __init__(self, 
        fill_resolution_low = None, 
        fill_resolution_high = None,
        out_suffix = ".filled.mtz",
        delete_tmp_files = True,
        ):

        self.fill_resolution_low = fill_resolution_low
        self.fill_resolution_high = fill_resolution_high
        self.out_suffix = out_suffix
        self.delete_tmp_files = delete_tmp_files

    def __call__(self,
        in_mtz, 
        ):
        """Complete the set of miller indices in an MTZ file"""

        logger.subheading('Filling missing reflections')

        out_mtz = in_mtz.with_suffix(self.out_suffix)

        tmp_mtz_1 = out_mtz.with_suffix('.step1-truncate.mtz')
        tmp_mtz_2 = out_mtz.with_suffix('.step2-uniquify.mtz')
        tmp_mtz_3 = out_mtz.with_suffix('.step3-remerged.mtz')

        fill_low, fill_high = self.get_fill_limits(in_mtz)

        # Stage 1 - truncate dataset
        self.truncate(
            in_mtz = in_mtz, 
            out_mtz = tmp_mtz_1,
            resolution_low = fill_low, 
            resolution_high = fill_high,
            )

        # Stage 2 - Uniqueify the file
        self.uniqueify(
            in_mtz = tmp_mtz_1, 
            out_mtz = tmp_mtz_2,
            )

        # Stage 3 - remerge the two files
        self.remerge_with_dummy_column(
            in_mtz_1 = in_mtz,
            in_mtz_2 = tmp_mtz_2,
            out_mtz = tmp_mtz_3,
            )

        # Stage 4 - remove the dummy column
        self.remove_dummy_column(
            in_mtz = tmp_mtz_3,
            out_mtz = out_mtz,
            )

        self.cleanup(
            tmp_files = [
                tmp_mtz_1, 
                tmp_mtz_2, 
                tmp_mtz_3, 
                ],
            )

        return out_mtz

    def get_fill_limits(self, 
        in_mtz,
        ):

        fill_low = self.fill_resolution_low
        fill_high = self.fill_resolution_high

        m = CrystallographicData.from_file(
            str(in_mtz)
            )

        if fill_low is None: 
            fill_low = m.crystal.resolution_low
        
        if fill_high is None: 
            fill_high = m.crystal.resolution_high

        if fill_low is None: 
            raise Exception('No fill_low supplied or extracted from mtz')

        if fill_high is None: 
            raise Exception('No fill_high supplied or extracted from mtz')

        return (fill_low, fill_high)

    def truncate(self, 
        in_mtz, 
        out_mtz, 
        resolution_low, 
        resolution_high,
        ):

        assert resolution_low is not None
        assert resolution_high is not None

        prog = Dispatcher('cad')

        prog.extend_args([
            'hklin1', str(in_mtz),
            'hklout', str(out_mtz),
        ])

        prog.extend_stdin([
            'monitor BRIEF',
            'labin file_number 1 ALL',
            'resolution file 1 {low} {high}'.format(
                low = resolution_low, 
                high = resolution_high,
                ),
        ])

        run_program(prog)

        if not out_mtz.exists():

            prog.write_output(
                str(out_mtz.with_suffix('.log'))
                )

            raise_missing(
                filepath = out_mtz, 
                result = prog.result,
                )

    def uniqueify(self, 
        in_mtz, 
        out_mtz,
        ):

        prog = Dispatcher('uniqueify')

        prog.extend_args([
            '-p', '0.05',
            str(in_mtz),
            str(out_mtz),
        ])

        run_program(prog)

        if not out_mtz.exists():

            prog.write_output(
                str(out_mtz.with_suffix('.log'))
                )

            raise_missing(
                filepath = out_mtz, 
                result = prog.result,
                )

        # Annoying uncontrollable log file! (in the cwd)
        possible_log = pl.Path(out_mtz.name).with_suffix('.log')
        if possible_log.exists():
            os.remove(str(possible_log))

    def remerge_with_dummy_column(self, in_mtz_1, in_mtz_2, out_mtz):

        prog = Dispatcher('cad')

        prog.extend_args([
            'hklin1', str(in_mtz_1),
            'hklin2', str(in_mtz_2),
            'hklout', str(out_mtz),
        ])

        prog.extend_stdin([
            'monitor BRIEF',
            'labin file_number 1 ALL',
            'labin file_number 2 E1=FreeR_flag',
            'labout file_number 2 E1=dummy',
        ])

        run_program(prog)

        if not out_mtz.exists():

            prog.write_output(
                str(out_mtz.with_suffix('.log'))
                )
            
            raise_missing(
                filepath = out_mtz, 
                result = prog.result,
                )

    def remove_dummy_column(self, in_mtz, out_mtz):

        prog = Dispatcher('mtzutils')

        prog.extend_args([
            'hklin1', str(in_mtz),
            'hklout', str(out_mtz),
        ])

        prog.extend_stdin([
            'HEADER BRIEF',
            'EXCLUDE 1 dummy',
            'ONEFILE',
            'END',
        ])

        run_program(prog)

        if not out_mtz.exists():

            prog.write_output(
                str(out_mtz.with_suffix('.log'))
                )
            
            raise_missing(
                filepath = out_mtz, 
                result = prog.result,
                )

    def cleanup(self, tmp_files):

        if (self.delete_tmp_files is True):
            for f in tmp_files:
                os.remove(str(f))


class TransferRFreeFlags(object):

    def __init__(self, 
        reference_mtz, 
        input_free_r_flag, 
        output_free_r_flag = None,
        out_suffix = ".free.mtz",
        delete_tmp_files = True,
        ):
        
        if output_free_r_flag is None:
            output_free_r_flag = input_free_r_flag

        self.reference_mtz = reference_mtz
        self.input_free_r_flag = input_free_r_flag 
        self.output_free_r_flag = output_free_r_flag
        self.out_suffix = out_suffix
        self.delete_tmp_files = delete_tmp_files

    def __call__(self,
        in_mtz, 
        ):
        """Copy R-free flags from reference mtz"""

        logger.subheading('Transferring and completing R-free flags')

        out_mtz = in_mtz.with_suffix(self.out_suffix)

        tmp_mtz = out_mtz.with_suffix('.step1-transfer.mtz')

        # Stage 1 - transfer R-free from reference
        self.transfer_rfree_flags(
            in_mtz = in_mtz,
            reference_mtz = self.reference_mtz,
            out_mtz = tmp_mtz,
            input_free_r_flag = self.input_free_r_flag,
            output_free_r_flag = self.output_free_r_flag,
            )

        # Stage 2 - populate missing R-free values
        self.fill_missing_flags(
            in_mtz = tmp_mtz,
            out_mtz = out_mtz,
            free_r_flag = self.output_free_r_flag,
            )

        if (self.delete_tmp_files is True):
            os.remove(str(tmp_mtz))

        return out_mtz

    def transfer_rfree_flags(self, 
        in_mtz,
        reference_mtz,
        out_mtz,
        input_free_r_flag,
        output_free_r_flag,
        ):

        prog = Dispatcher('cad')

        prog.extend_args([
            'hklin1', str(in_mtz),
            'hklin2', str(reference_mtz),
            'hklout', str(out_mtz),
        ])

        prog.extend_stdin([
            'labin file_number 1 ALL',
            'labin file_number 2 E1={}'.format(
                input_free_r_flag
                ),
            'labout file_number 2 E1={}'.format(
                output_free_r_flag
                ),
            'END',
        ])

        run_program(prog)

        if not out_mtz.exists():

            prog.write_output(
                str(out_mtz.with_suffix('.log'))
                )
            
            raise_missing(
                filepath = out_mtz, 
                result = prog.result,
                )

    def fill_missing_flags(self, 
        in_mtz, 
        out_mtz, 
        free_r_flag,
        ):

        prog = Dispatcher('freerflag')

        prog.extend_args([
            'hklin',  str(in_mtz),
            'hklout', str(out_mtz),
        ])

        prog.extend_stdin([
            'COMPLETE FREE={}'.format(free_r_flag),
            'END',
        ])

        run_program(prog)

        if not out_mtz.exists():

            prog.write_output(
                str(out_mtz.with_suffix('.log'))
                )
            
            raise_missing(
                filepath = out_mtz, 
                result = prog.result,
                )


class PrepareMTZTask(object):

    def __init__(self, 
        remove_columns = None,
        reindex_to_reference = None,
        fill_missing_reflections = None,
        transfer_rfree_flags = None,
        output_directory = None,
        path_labeller = None,
        delete_tmp_files = True,
        ):

        self.remove_columns = (
            remove_columns
            )
        
        self.reindex_to_reference = (
            reindex_to_reference
            )

        self.fill_missing_reflections = (
            fill_missing_reflections
            )

        self.transfer_rfree_flags = (
            transfer_rfree_flags
            )

        self.output_directory = (
            output_directory
            )

        self.path_labeller = (
            path_labeller 
            if 
            (path_labeller is not None)
            else 
            PathLabeller()
            )

        self.delete_tmp_files = (
            delete_tmp_files
            )

    def __call__(self, in_mtz):
        """Prepare mtz by various standard sanitisations"""

        in_mtz = pl.Path(in_mtz)

        label = self.path_labeller(
            str(in_mtz)
            )

        if self.output_directory is not None:

            out_dir = (
                self.output_directory / label
                )

            if not out_dir.is_dir():
                out_dir.mkdir(parents=True)

            # Make a copy of the file in the new directory
            in_mtz = self.copy_to_directory(
                filepath = in_mtz,
                directory = out_dir,
                )

        else: 

            # All done in place
            out_dir = None 

        file_handler, warning_handler = self.start_logging(
            logpath = in_mtz.with_suffix('.prepared.log'),
            )

        logger.subheading(
            'Preparing MTZ for {}'.format(
                label,
                )
            )

        try: 

            out_mtz = self.run(
                in_mtz = in_mtz, 
                )

        except Exception as e: 

            # Ensure usable information is logged in the local logfile
            logger.warning(traceback.format_exc())

            self.stop_logging(
                handlers = [file_handler, warning_handler],
                )

            return TaskReturn(
                output = ModelDataInputOutput(
                    input_mtz = in_mtz,
                    ),
                status = TaskReturnStatus(
                    success = False,
                    errors = [e],
                    warnings = list(warning_handler.list()),
                    ),
                )

        self.stop_logging(
            handlers = [file_handler, warning_handler],
            )

        return TaskReturn(
            output = ModelDataInputOutput(
                input_mtz = in_mtz,
                output_mtz = out_mtz,
                ),
            status = TaskReturnStatus(
                success = True,
                errors = None,
                warnings = list(warning_handler.list()),
                ),
            )

    def run(self, in_mtz):

        # Save the original mtz as variable
        orig_mtz = in_mtz

        # List of intermediate files to be deleted at the end
        intermediate_files = []

        if self.remove_columns is not None: 

            if str(orig_mtz) != str(in_mtz):
                intermediate_files.append(in_mtz)

            in_mtz = self.remove_columns(
                in_mtz = in_mtz,
                )

        if self.reindex_to_reference is not None:
            
            if str(orig_mtz) != str(in_mtz):
                intermediate_files.append(in_mtz)

            in_mtz = self.reindex_to_reference(
                in_mtz = in_mtz,
            )

        if self.fill_missing_reflections is not None:
            
            if str(orig_mtz) != str(in_mtz):
                intermediate_files.append(in_mtz)

            in_mtz = self.fill_missing_reflections(
                in_mtz = in_mtz,
                )

        if self.transfer_rfree_flags is not None: 

            if str(orig_mtz) != str(in_mtz):
                intermediate_files.append(in_mtz)

            in_mtz = self.transfer_rfree_flags(
                in_mtz = in_mtz,
                )

        ###

        if str(orig_mtz) != str(in_mtz):

            out_mtz = orig_mtz.with_suffix('.prepared.mtz')

            shutil.copy(
                str(in_mtz), 
                str(out_mtz),
                )
            
            intermediate_files.append(in_mtz)

        else: 

            # If the input file is unchanged -- return input file
            out_mtz = orig_mtz

        # Tidy up
        if self.delete_tmp_files:
            for f in intermediate_files:
                # Do not allow to delete input
                if str(f) == str(orig_mtz):
                    continue
                else:
                    os.remove(str(f))

        return out_mtz

    def start_logging(self, logpath):

        logger = lg.getLogger('giant')

        f_handler = lg.lg.FileHandler(
            filename = str(logpath),
            )

        logger.addHandler(f_handler)

        w_handler = lg.add_warning_handler(
            logger = logger,
            )

        return f_handler, w_handler

    def stop_logging(self, handlers):

        logger = lg.getLogger('giant')

        for h in handlers:
            logger.removeHandler(h)

    def copy_to_directory(self, filepath, directory):

        out_path = (
            directory / filepath.name
            )

        shutil.copy(
            str(filepath), 
            str(out_path),
            )

        return out_path
