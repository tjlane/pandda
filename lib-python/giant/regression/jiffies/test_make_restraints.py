import giant.logs as lg
logger = lg.getLogger(__name__)

import os, shutil

from pytest import mark

from iotbx.pdb.atom_selection import AtomSelectionError

from giant.regression.jiffies import run_jiffy, check_files_exist

def test_make_restraints(tmpdir):

    import giant.jiffies.make_restraints as make_r

    from giant.regression.jiffies.multi_conf_test_data import \
            LIGCODE, MERGED_PDB_STR

    in_dir = tmpdir

    # Write input string to PDB file
    input_pdb = os.path.join(str(in_dir), 'input.pdb')
    with open(input_pdb, 'w') as fh:
        fh.write(MERGED_PDB_STR)

    dp = make_r.master_phil.extract()

    output_log = os.path.join(str(in_dir), 'make_restraints.log')
    output_restraints_phenix = os.path.join(str(in_dir), 'phenix.params')
    output_restraints_refmac = os.path.join(str(in_dir), 'refmac.params')

    make_args = [
        input_pdb,
        'output.log={}'.format(output_log),
        'occupancy.resname={}'.format(LIGCODE),
        "output.phenix={}".format(output_restraints_phenix),
        "output.refmac={}".format(output_restraints_refmac),
    ]

    expected_output_files = [
        output_log,
        output_restraints_phenix,
        output_restraints_refmac,
    ]

    run_jiffy(
        args = make_args,
        module = make_r,
    )

    check_files_exist(expected_output_files)

    with open(output_restraints_phenix) as fh:
        assert fh.read().strip() == "\n".join([
            "refinement.refine.occupancies {",
            "    constrained_group {",
            "        selection = (chain 'A' and resseq    1 and icode ' ' and resname 'HOH' and altid 'A') or \\",
            "                    (chain 'A' and resseq    2 and icode ' ' and resname 'HOH' and altid 'A') or \\",
            "                    (chain 'A' and resseq    4 and icode ' ' and resname 'HOH' and altid 'A') or \\",
            "                    (chain 'B' and resseq   53 and icode ' ' and resname 'ARG' and altid 'A') or \\",
            "                    (chain 'B' and resseq   54 and icode ' ' and resname 'ILE' and altid 'A') or \\",
            "                    (chain 'B' and resseq  116 and icode ' ' and resname 'PHE' and altid 'A')",
            "        selection = (chain 'A' and resseq    4 and icode ' ' and resname 'HOH' and altid 'B') or \\",
            "                    (chain 'B' and resseq   53 and icode ' ' and resname 'ARG' and altid 'B') or \\",
            "                    (chain 'B' and resseq   54 and icode ' ' and resname 'ILE' and altid 'B') or \\",
            "                    (chain 'B' and resseq  116 and icode ' ' and resname 'PHE' and altid 'B') or \\",
            "                    (chain 'B' and resseq 1145 and icode ' ' and resname 'PW3' and altid 'B')",
            "    }",
            "}",
        ]).strip()

    with open(output_restraints_refmac) as fh:
        assert fh.read().strip() == "\n".join([
            "occupancy group id 1 chain A resi    1 alte A",
            "occupancy group id 1 chain A resi    2 alte A",
            "occupancy group id 1 chain A resi    4 alte A",
            "occupancy group id 1 chain B resi   53 alte A",
            "occupancy group id 1 chain B resi   54 alte A",
            "occupancy group id 1 chain B resi  116 alte A",
            "occupancy group id 2 chain A resi    4 alte B",
            "occupancy group id 2 chain B resi   53 alte B",
            "occupancy group id 2 chain B resi   54 alte B",
            "occupancy group id 2 chain B resi  116 alte B",
            "occupancy group id 2 chain B resi 1145 alte B",
            "occupancy group alts complete 1 2",
            "occupancy refine",
        ]).strip()

