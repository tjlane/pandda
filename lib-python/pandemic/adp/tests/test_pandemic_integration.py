from pytest import approx, raises, mark
from libtbx import adopt_init_args
from functools import partial
import itertools

from pandemic.adp import run_pandemic_adp

from giant.paths import filename, foldername, not_available

# Programs required for tests that are not available
MISSING_PROGRAMS = tuple(not_available(['pymol']))

def fmt_option(option_name, option_value):
    return "{!s}={!s}".format(option_name, option_value)

def fmt_out_dir(output_directory):
    return "output.out_dir={}".format(output_directory)

def make_options(option_name, option_choices):
    return [fmt_option(option_name, c) for c in option_choices]

def combine_options(options_list):
    return itertools.product(*options_list)

def setup_input_structures(output_directory, number_of_structures, include_cryst_lines=True):

    from pandemic.resources.structure_test_snippets import \
        pdb_test_structure_atoms, \
        pdb_test_structure_cryst_lines
    assert len(pdb_test_structure_atoms) == len(pdb_test_structure_cryst_lines)

    structure_strings = itertools.cycle(list(zip(pdb_test_structure_cryst_lines, pdb_test_structure_atoms)))

    structure_directory = output_directory.mkdir('structures')

    output_filenames = []

    for i in range(number_of_structures):

        s_directory = structure_directory.mkdir('s{}'.format(i+1))
        s_filename = s_directory / 'structure{}.pdb'.format(i+1)

        cryst_string, atoms_string = next(structure_strings)

        if include_cryst_lines is True:
            s_filename.write(cryst_string+'\n', mode='a')
        s_filename.write(atoms_string+'\n', mode='a')

        output_filenames.append(str(s_filename))

    return output_filenames


class _CheckPandemicOutputFiles(object):

    chains = None

    def __init__(self,
        tls_level_names,
        labelling,
        pymol_images,
        ):

        # Overwrite if pymol not installed
        if 'pymol' in MISSING_PROGRAMS:
            pymol_images = False

        # How many levels are expected?
        n_tls_levels = len(tls_level_names)
        n_all_levels = len(tls_level_names) + 1
        # Numbers of each level
        tls_level_numbers = tuple([i+1 for i in range(n_tls_levels)])
        adp_level_number  = n_all_levels
        all_level_numbers = tuple([i+1 for i in range(n_all_levels)])

        adopt_init_args(self, locals())

        self.init_expected_files()

    def get_structure_labels(self, filenames):

        if self.labelling == 'filename':
            return [filename(f) for f in filenames]
        elif self.labelling == 'foldername':
            return [foldername(f) for f in filenames]
        else:
            raise Exception('bad test')

    def check_all(self,
        output_directory,
        input_structures,
        ):

        self.check_misc(
            output_directory=output_directory,
            )

        self.check_level_partitions(
            output_directory=output_directory,
            )

        self.check_hierarchical_model(
            output_directory=output_directory,
            )

        self.check_optimisation(
            output_directory=output_directory,
            )

        self.check_output_structures(
            output_directory=output_directory,
            input_structures=input_structures,
            )

        self.check_analysis(
            output_directory=output_directory,
            )

        self.check_pymol_images(
            output_directory=output_directory,
            )

    def check_misc(self, output_directory):

        for f in self.misc_files:
            self.check_file(output_directory/f)

    def check_level_partitions(self, output_directory):

        for f in self.level_partition_files:
            self.check_file(output_directory/f)

    def check_hierarchical_model(self, output_directory):

        for f in self.hierarchical_model_files:
            self.check_file(output_directory/f)

    def check_optimisation(self, output_directory):

        for f in self.optimisation_files:
            self.check_file(output_directory/f)

    def check_output_structures(self, output_directory, input_structures):

        structure_labels = self.get_structure_labels(input_structures)

        for s in structure_labels:
            for f in self.structures_files:
                self.check_file(output_directory/f.format(label=s))

    def check_analysis(self, output_directory):

        for f in self.analysis_files:
            self.check_file(output_directory/f)

    def check_pymol_images(self, output_directory):

        for f in self.pymol_files:
            self.check_file(output_directory/f, expected_exists=self.pymol_images)

    def check_file(self, filename, expected_exists=True):
        if not (filename.exists() is expected_exists):
            raise Exception('File does/does not exist: {}'.format(str(filename)))


class CheckPandemicOutputFiles(_CheckPandemicOutputFiles):

    chains = ['A']

    def init_expected_files(self):

        self.misc_files = [
            "pandemic.log",
            "params-input.eff",
            "params-running.eff",
            "results.html",
            "output_data.csv",
            "starting_model.json",
            "optimised_model.json",
        ]

        ###################################################

        self.level_partition_files = [
            "level_partitions/pymol_script.py",
        ]

        self.hierarchical_model_files = [
            "hierarchical_model/pymol_script.py",
            "hierarchical_model/structures/average_uijs_input.pdb",
            "hierarchical_model/structures/average_uijs_output.pdb",
        ]

        for chain in self.chains:
            for f in [
                # currently broken
                #"hierarchical_model/anisotropy-all-levels-chain_{chain}.png",
                "hierarchical_model/profile-all-levels-chain_{chain}.png",
                ]:
                self.hierarchical_model_files.append(f.format(chain=chain))

        for leveln in self.tls_level_numbers:
            for f in [
                "hierarchical_model/tables/tls_amplitudes_level_{leveln:04d}.csv",
                "hierarchical_model/tables/tls_matrices_level_{leveln:04d}.csv",
                "hierarchical_model/tables/tls_origins_level_{leveln:04d}.csv",
                "hierarchical_model/structures/level_{leveln}-tls_mode_1.pdb",
                ]:
                self.hierarchical_model_files.append(f.format(leveln=leveln))
        for leveln in self.all_level_numbers:
            for f in [
                "hierarchical_model/structures/level_{leveln}-all.pdb",
                ]:
                self.hierarchical_model_files.append(f.format(leveln=leveln))

        for chain, leveln in itertools.product(self.chains, self.all_level_numbers):
            for f in [
                # currently broken
                #"hierarchical_model/anisotropy-level_{leveln}-chain_{chain}.png",
                #"hierarchical_model/anisotropy-level_{leveln}-chain_{chain}.png",
                #"hierarchical_model/anisotropy-level_{leveln}-chain_{chain}.png",
                "hierarchical_model/profile-level_{leveln}-chain_{chain}.png",
                "hierarchical_model/profile-level_{leveln}-chain_{chain}.png",
                "hierarchical_model/profile-level_{leveln}-chain_{chain}.png",
                ]:
                self.hierarchical_model_files.append(f.format(leveln=leveln, chain=chain))

        ###################################################

        self.optimisation_files = [
            "optimisation/level_amplitudes_weights-sum_of_amplitudes.png",
            "optimisation/level_amplitudes_weights-sum_of_amplitudes_squared.png",
            "optimisation/level_amplitudes_weights-sum_of_squared_amplitudes.png",
            "optimisation/tracking_amplitudes.png",
            "optimisation/tracking_convergence.png",
            "optimisation/tracking_snapshots.png",
            "optimisation/tracking_rmsds.png",
            "optimisation/tracking_echt.csv",
            "optimisation/tracking_levels.csv",
            "optimisation/tracking_rmsds.csv",
        ]
        for chain in self.chains:
            for f in [
                "optimisation/tracking_atoms-chain_{chain}.png",
                ]:
                self.optimisation_files.append(f.format(chain=chain))

        ###################################################

        self.structures_files = [
            "structures/{label}/pymol_script.py",
            "structures/{label}/{label}.all-levels.pdb",
            "structures/{label}/{label}.all-tls-levels.pdb",
            "structures/{label}/{label}.all.pdb",
            "structures/{label}/{label}.atomic-level.pdb",
            "structures/{label}/{label}.input.pdb",
        ]

        for leveln in self.tls_level_numbers:
            for f in [
                "structures/{{label}}/{{label}}.tls-level-{leveln:04d}.pdb",
                ]:
                self.structures_files.append(f.format(leveln=leveln))

        for leveli, levelj in itertools.product(self.tls_level_numbers, self.tls_level_numbers):
            if leveli >= levelj: continue
            for f in [
                "structures/{{label}}/{{label}}.tls-level-{leveli:04d}-to-{levelj:04d}.pdb",
                ]:
                self.structures_files.append(f.format(leveli=leveli, levelj=levelj))

        ###################################################

        self.analysis_files = [
            "analysis/hierarchy_groups/input_hierarchy.eff",
            "analysis/hierarchy_groups/output_hierarchy.eff",
            "analysis/residuals/residual_vs_bfactor.png",
            # specific ones to this test set -- refactor in future?
            "analysis/residuals/residual_by_atom_type_1-30.png",
            "analysis/residuals/residual_by_atom_type_31-60.png",
            "analysis/residuals/residual_by_atom_type_61-90.png",
            "analysis/residuals/residual_by_atom_type_91-120.png",
            "analysis/residuals/residual_by_atom_type_121-150.png",
            "analysis/tls_group_clustering/tls_clustering_level_2-new_levels.png",
        ]

        for chain in self.chains:
            for f in [
                "analysis/hierarchy_groups/output-level-partitions-chain-{chain}.png",
                "analysis/residuals/residual_by_residue-chain_{chain}.png",
                ]:
                self.analysis_files.append(f.format(chain=chain))

        for chain, leveln in itertools.product(self.chains, self.all_level_numbers):
            for f in [
                "analysis/residuals/fitting_residual_correlations_level_{leveln}-chain_{chain}.png",
                ]:
                self.analysis_files.append(f.format(chain=chain, leveln=leveln))

        ###################################################

        self.pymol_files = [

        ]


class TestPandemicIntegrationSingleStructure(object):

    short_run_options = [
        "max_macro_cycles=1",
        "max_micro_cycles=1",
        "auto_levels=chain+ss",
        "atomic_adp_level=False",
        "dpi=50",
        "embed_images=False",
        "pymol=none",
        ]

    dry_run = [
        "dry_run=True",
        ]

    def run_simple_options_lists(self,
        tmpdir,
        options_lists,
        cryst_line=True,
        dry_run=False,
        ):

        input_structures = setup_input_structures(tmpdir, 1, include_cryst_lines=cryst_line)

        output_directories = []

        for options in options_lists:
            run_directory = tmpdir.make_numbered_dir(prefix='test-', rootdir=tmpdir) / 'pandemic-adp'
            output_directories.append(run_directory)
            run_options = [
                fmt_out_dir(run_directory),
                ] + \
                self.short_run_options + \
                (self.dry_run if (dry_run is True) else []) + \
                input_structures + \
                options
            run_pandemic_adp(run_options)

        return input_structures, output_directories

    def run_simple_option_check(self,
        tmpdir,
        option_name,
        option_values,
        cryst_line=True,
        dry_run=False,
        ):

        options_lists = [[fmt_option(option_name, v)] for v in option_values]

        return self.run_simple_options_lists(
            tmpdir=tmpdir,
            options_lists=options_lists,
            cryst_line=cryst_line,
            dry_run=dry_run,
            )

    @mark.slow
    def test_no_args(self, tmpdir):

        # Need to run with something so choose html (default on anyway)
        structures, directories = self.run_simple_option_check(
            tmpdir=tmpdir,
            option_name='make_html',
            option_values=[True],
            cryst_line=True,
            dry_run=False,
            )

        checker = CheckPandemicOutputFiles(
            tls_level_names=['chain','ss'],
            labelling='filename',
            pymol_images=True,
            )

        checker.check_all(
            output_directory=directories[0],
            input_structures=structures,
            )

    def test_file_labelling(self, tmpdir):

        with raises(SystemExit) as e:
            self.run_simple_option_check(
                tmpdir=tmpdir,
                option_name='input.labelling',
                option_values=['filename','foldername'],
                cryst_line=True,
                dry_run=True,
                )
        assert str(e.value) == 'Program exited normally'

    @mark.slow
    def test_model_type(self, tmpdir):

        self.run_simple_option_check(
            tmpdir=tmpdir.mkdir('check1'),
            option_name='model_type',
            option_values=['crystallographic'],
            cryst_line=True,
            dry_run=False,
            )

        self.run_simple_option_check(
            tmpdir=tmpdir.mkdir('check2'),
            option_name='model_type',
            option_values=['noncrystallographic'],
            cryst_line=False,
            dry_run=False,
            )

    def test_model_type_failures(self, tmpdir):

        # Should not error (automatically identifies non-crystallographic models)
        with raises(SystemExit) as e:
            self.run_simple_option_check(
                tmpdir=tmpdir.mkdir('check1'),
                option_name='model_type',
                option_values=['crystallographic'],
                cryst_line=False,
                dry_run=True,
                )
        assert str(e.value) == 'Program exited normally'
        #msg = "Failed to load pdb file -- have you selected the correct model_type? Current model_type is crystallographic."
        #assert str(e.value).startswith(msg)

        # Shouldn't care.
        with raises(SystemExit) as e:
            self.run_simple_option_check(
                tmpdir=tmpdir.mkdir('check2'),
                option_name='model_type',
                option_values=['noncrystallographic'],
                cryst_line=True,
                dry_run=True,
                )
        assert str(e.value) == 'Program exited normally'

    @mark.slow
    def test_pymol(self, tmpdir):

        # This should run regardless of whether pymol is on the system
        self.run_simple_option_check(
            tmpdir=tmpdir.mkdir('check1'),
            option_name='pymol',
            option_values=['none'],
            cryst_line=True,
            dry_run=False,
            )
        # TODO: CHECK FILES DON'T EXIST

        # Create prepped function that can be used regardless of whether pymol is installed
        run_with_pymol = partial(
            self.run_simple_option_check,
            tmpdir=tmpdir.mkdir('check2'),
            option_name='pymol',
            option_values=['none'],
            cryst_line=True,
            dry_run=False,
            )

        if 'pymol' in MISSING_PROGRAMS:
            with raises(Exception) as e:
                run_with_pymol()
            assert """pymol is required when output.images.pymol is not set to "none".""" in str(e.value)
            assert """The program is not currently available in any of the PATH directories.""" in str(e.value)
            assert """[ It must be available as an executable script/link -- not as an alias.  ]""" in str(e.value)
            assert """To turn off pymol add output.images.pymol=none to the command line options.""" in str(e.value)
        else:
            structures, directories = run_with_pymol()

            # Check existence of output files
            checker = CheckPandemicOutputFiles(
                tls_level_names=['chain','ss'],
                labelling='filename',
                pymol_images=True,
                )

            for d in directories:
                checker.check_pymol_images(
                    output_directory=d,
                    )

