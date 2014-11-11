

apo_struct = '/home/npearce/LIGANDTEST/x529/1-apo/apo-refmac.pdb'
apo_data = '/home/npearce/LIGANDTEST/x529/1-apo/apo-refmac.mtz'

ligand_pdb = '/home/npearce/LIGANDTEST/x529/2-ligand/ligand-grade.pdb'
ligand_cif = '/home/npearce/LIGANDTEST/x529/2-ligand/ligand-grade.cif'

# ====

print('=====>\nvariables set!\n=====>')

imol_map = make_and_draw_map(apo_data, "FWT", "PHWT", "", 0, 0)
imol_protein = read_pdb(apo_struct)

# don't do this
###set_rotation_centre(54, 10, 20)

ligand_mol = handle_read_draw_molecule_with_recentre(ligand_pdb, 0)
read_cif_dictionary(ligand_cif)

print('=====>\nimported!\n=====>')

set_find_ligand_here_cluster(1)
set_ligand_search_protein_molecule(imol_protein)
set_ligand_search_map_molecule(imol_map)
add_ligand_search_wiggly_ligand_molecule(ligand_mol)

print('=====>\nligandfitting ready to go!\n=====>')

results = execute_ligand_search_py()

print("=============================================")
print(results)
print("=============================================")

for result in results:
    refine_residues(result, [['A', 1, '']])
    accept_regularizement()
    ligand_file_name = "fitted-ligand-" + str(result) + ".pdb"
    write_pdb_file(result, ligand_file_name)




#comp_id = "3GP"
#
#imol_map = make_and_draw_map("refmac.mtz", "FWT", "PHWT", "", 0, 0)
#imol_protein = read_pdb("refmac.pdb")
#
## don't do this
#set_rotation_centre(54, 10, 20)
#
#
#pdb_file_name  = comp_id + ".pdb"
#dict_file_name = comp_id + ".cif"
#
#ligand_mol = handle_read_draw_molecule_with_recentre(pdb_file_name, 0)
#read_cif_dictionary(dict_file_name)
#
#set_find_ligand_here_cluster_flag(1)
#set_ligand_search_protein_molecule(imol_protein)
#set_ligand_search_map_molecule(imol_map)
#add_ligand_search_wiggly_ligand_molecule(ligand_mol)
#
#results = execute_ligand_search_py()
#
#print "============================================="
#print results
#print "============================================="
#
#for result in results:
#    refine_residues(result, [['A', 1, '']])
#    accept_regularizement()
#    ligand_file_name = "fitted-ligand-" + str(result) + ".pdb"
#    write_pdb_file(result, ligand_file_name)
#
