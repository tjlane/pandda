import sys

__name__ = "pandda_test"

execfile("/home/npearce/bin/PANDDAs/bin/pandda.analyse")

input_args = [
"testing=True",
"cpus=8",
"verbose=True",
#
#"data_dirs='/home/npearce/Analyses/FragmentSoaks/BAZ2BA-Feb13-Soak/ProcessedFragmentSoak/BAZ2BA-x52*/1-apo'",
#"pdb_style='apo-BAZ2BA-x52*-refmac.pdb'",
#
"data_dirs='/home/npearce/Analyses/FragmentSoaks/BAZ2BA-Feb13-Soak/ProcessedFragmentSoak/BAZ2BA-*/1-apo'",
"pdb_style='apo-BAZ2BA-*-refmac.pdb'",
#
"lig_style='../2-ligand/*.cif'",
#
"outdir=/home/npearce/Analyses/FragmentSoaks/BAZ2BA-Feb13-Soak/pandda-2.5-local-test",
"alignment.method=local",
"contour_level=2.5",
"resolution_factor=0.33",
#
"maps.ampl_label=FWT",
"maps.phas_label=PHWT"
]

pandda=pandda_main(input_args)
