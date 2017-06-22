


original = """
HETATM 2053  C01 PW3 B1145     -28.613   3.252   8.415  0.87 29.90           C
HETATM 2054  C02 PW3 B1145     -28.348   3.089   6.927  0.87 29.23           C
HETATM 2055  N03 PW3 B1145     -29.447   2.886   6.061  0.87 38.58           N
HETATM 2056  C04 PW3 B1145     -30.772   2.819   6.612  0.87 41.88           C
HETATM 2057  C05 PW3 B1145     -31.838   2.590   5.636  0.87 42.90           C
HETATM 2058  N06 PW3 B1145     -31.685   3.422   4.449  0.87 35.35           N
HETATM 2059  C07 PW3 B1145     -30.380   3.510   3.946  0.87 38.37           C
HETATM 2060  C08 PW3 B1145     -29.295   2.704   4.713  0.87 35.02           C
HETATM 2061  O09 PW3 B1145     -30.116   4.169   2.970  0.87 39.95           O
HETATM 2062  O10 PW3 B1145     -27.236   3.110   6.465  0.87 25.71           O
"""

unsorted = """
HETATM 2053  C01 PW3 B1145     -28.613   3.252   8.415  0.87 29.90           C
HETATM 2054  C05 PW3 B1145     -31.838   2.590   5.636  0.87 42.90           C
HETATM 2055  N06 PW3 B1145     -31.685   3.422   4.449  0.87 35.35           N
HETATM 2056  C07 PW3 B1145     -30.380   3.510   3.946  0.87 38.37           C
HETATM 2057  C08 PW3 B1145     -29.295   2.704   4.713  0.87 35.02           C
HETATM 2058  C02 PW3 B1145     -28.348   3.089   6.927  0.87 29.23           C
HETATM 2059  N03 PW3 B1145     -29.447   2.886   6.061  0.87 38.58           N
HETATM 2060  C04 PW3 B1145     -30.772   2.819   6.612  0.87 41.88           C
HETATM 2061  O09 PW3 B1145     -30.116   4.169   2.970  0.87 39.95           O
HETATM 2062  O10 PW3 B1145     -27.236   3.110   6.465  0.87 25.71           O
"""


h_orig = iotbx.pdb.hierarchy.input(pdb_string=original).hierarchy
h_mixd = iotbx.pdb.hierarchy.input(pdb_string=unsorted).hierarchy

h_sort = h_mixd.deep_copy()
h_sort.sort_atoms_in_place()

