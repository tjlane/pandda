from pytest import approx, raises, mark

from libtbx import group_args
import iotbx.pdb

from pandemic.resources.structure_test_snippets import pdb_test_structure_atoms

def test_GenerateLevelSelectionsTask_auto_levels_1():

    from pandemic.adp.hierarchy.level_selections import GenerateLevelSelectionsTask

    hierarchy = iotbx.pdb.hierarchy.input(pdb_string=pdb_test_structure_atoms[0]).hierarchy

    for input_levels, level_labels, correct_groups in [
        [
            ['chain', 'ss', 'residue', 'backbone_sidechain'],
            ['chain', 'sec. struct.', 'residue', 'backbone/sidechain'],
            {
                'chain' : ["chain 'A'"],
                'sec. struct.' : [
                    "chain 'A' and resid 1900  through 1905",
                    "chain 'A' and resid 1906  through 1909",
                    "chain 'A' and resid 1910  through 1920",
                    "chain 'A' and resid 1921  through 1924",
                    "chain 'A' and resid 1925  through 1944",
                    "chain 'A' and resid 1945  through 1950",
                    ],
                'residue' : [
                    "chain 'A' and resseq 1900 and icode ' '",
                    "chain 'A' and resseq 1901 and icode ' '",
                    "chain 'A' and resseq 1902 and icode ' '",
                    "chain 'A' and resseq 1903 and icode ' '",
                    "chain 'A' and resseq 1904 and icode ' '",
                    "chain 'A' and resseq 1905 and icode ' '",
                    "chain 'A' and resseq 1906 and icode ' '",
                    "chain 'A' and resseq 1907 and icode ' '",
                    "chain 'A' and resseq 1908 and icode ' '",
                    "chain 'A' and resseq 1909 and icode ' '",
                    "chain 'A' and resseq 1910 and icode ' '",
                    "chain 'A' and resseq 1911 and icode ' '",
                    "chain 'A' and resseq 1912 and icode ' '",
                    "chain 'A' and resseq 1913 and icode ' '",
                    "chain 'A' and resseq 1914 and icode ' '",
                    "chain 'A' and resseq 1915 and icode ' '",
                    "chain 'A' and resseq 1916 and icode ' '",
                    "chain 'A' and resseq 1917 and icode ' '",
                    "chain 'A' and resseq 1918 and icode ' '",
                    "chain 'A' and resseq 1919 and icode ' '",
                    "chain 'A' and resseq 1920 and icode ' '",
                    "chain 'A' and resseq 1921 and icode ' '",
                    "chain 'A' and resseq 1922 and icode ' '",
                    "chain 'A' and resseq 1923 and icode ' '",
                    "chain 'A' and resseq 1924 and icode ' '",
                    "chain 'A' and resseq 1925 and icode ' '",
                    "chain 'A' and resseq 1926 and icode ' '",
                    "chain 'A' and resseq 1927 and icode ' '",
                    "chain 'A' and resseq 1928 and icode ' '",
                    "chain 'A' and resseq 1929 and icode ' '",
                    "chain 'A' and resseq 1930 and icode ' '",
                    "chain 'A' and resseq 1931 and icode ' '",
                    "chain 'A' and resseq 1932 and icode ' '",
                    "chain 'A' and resseq 1933 and icode ' '",
                    "chain 'A' and resseq 1934 and icode ' '",
                    "chain 'A' and resseq 1935 and icode ' '",
                    "chain 'A' and resseq 1936 and icode ' '",
                    "chain 'A' and resseq 1937 and icode ' '",
                    "chain 'A' and resseq 1938 and icode ' '",
                    "chain 'A' and resseq 1939 and icode ' '",
                    "chain 'A' and resseq 1940 and icode ' '",
                    "chain 'A' and resseq 1941 and icode ' '",
                    "chain 'A' and resseq 1942 and icode ' '",
                    "chain 'A' and resseq 1943 and icode ' '",
                    "chain 'A' and resseq 1944 and icode ' '",
                    "chain 'A' and resseq 1945 and icode ' '",
                    "chain 'A' and resseq 1946 and icode ' '",
                    "chain 'A' and resseq 1947 and icode ' '",
                    "chain 'A' and resseq 1948 and icode ' '",
                    "chain 'A' and resseq 1949 and icode ' '",
                    "chain 'A' and resseq 1950 and icode ' '",
                    ],
                'backbone/sidechain' : [
                    "chain 'A' and resseq 1901 and icode ' ' and resname 'TYR' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1902 and icode ' ' and resname 'LYS' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1903 and icode ' ' and resname 'LYS' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1904 and icode ' ' and resname 'VAL' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1905 and icode ' ' and resname 'ILE' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1906 and icode ' ' and resname 'LYS' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1907 and icode ' ' and resname 'LYS' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1909 and icode ' ' and resname 'MET' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1910 and icode ' ' and resname 'ASP' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1911 and icode ' ' and resname 'PHE' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1912 and icode ' ' and resname 'SER' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1913 and icode ' ' and resname 'THR' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1914 and icode ' ' and resname 'ILE' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1915 and icode ' ' and resname 'ARG' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1916 and icode ' ' and resname 'GLU' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1917 and icode ' ' and resname 'LYS' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1918 and icode ' ' and resname 'LEU' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1919 and icode ' ' and resname 'SER' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1920 and icode ' ' and resname 'SER' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1922 and icode ' ' and resname 'GLN' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1923 and icode ' ' and resname 'TYR' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1925 and icode ' ' and resname 'ASN' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1926 and icode ' ' and resname 'LEU' and altid 'A' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1926 and icode ' ' and resname 'LEU' and altid 'B' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1927 and icode ' ' and resname 'GLU' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1928 and icode ' ' and resname 'THR' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1929 and icode ' ' and resname 'PHE' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1931 and icode ' ' and resname 'LEU' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1932 and icode ' ' and resname 'ASP' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1933 and icode ' ' and resname 'VAL' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1934 and icode ' ' and resname 'ARG' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1935 and icode ' ' and resname 'LEU' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1936 and icode ' ' and resname 'VAL' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1937 and icode ' ' and resname 'PHE' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1938 and icode ' ' and resname 'ASP' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1939 and icode ' ' and resname 'ASN' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1940 and icode ' ' and resname 'CYS' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1941 and icode ' ' and resname 'GLU' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1942 and icode ' ' and resname 'THR' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1943 and icode ' ' and resname 'PHE' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1944 and icode ' ' and resname 'ASN' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1945 and icode ' ' and resname 'GLU' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1945 and icode ' ' and resname 'GLU' and altid 'A' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1945 and icode ' ' and resname 'GLU' and altid 'B' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1946 and icode ' ' and resname 'ASP' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1947 and icode ' ' and resname 'ASP' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1948 and icode ' ' and resname 'SER' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1949 and icode ' ' and resname 'ASP' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1950 and icode ' ' and resname 'ILE' and altid ' ' and (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1901 and icode ' ' and resname 'TYR' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1902 and icode ' ' and resname 'LYS' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1903 and icode ' ' and resname 'LYS' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1904 and icode ' ' and resname 'VAL' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1905 and icode ' ' and resname 'ILE' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1906 and icode ' ' and resname 'LYS' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1907 and icode ' ' and resname 'LYS' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1909 and icode ' ' and resname 'MET' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1910 and icode ' ' and resname 'ASP' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1911 and icode ' ' and resname 'PHE' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1912 and icode ' ' and resname 'SER' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1913 and icode ' ' and resname 'THR' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1914 and icode ' ' and resname 'ILE' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1915 and icode ' ' and resname 'ARG' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1916 and icode ' ' and resname 'GLU' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1917 and icode ' ' and resname 'LYS' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1918 and icode ' ' and resname 'LEU' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1919 and icode ' ' and resname 'SER' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1920 and icode ' ' and resname 'SER' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1922 and icode ' ' and resname 'GLN' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1923 and icode ' ' and resname 'TYR' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1925 and icode ' ' and resname 'ASN' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1926 and icode ' ' and resname 'LEU' and altid 'A' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1926 and icode ' ' and resname 'LEU' and altid 'B' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1927 and icode ' ' and resname 'GLU' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1928 and icode ' ' and resname 'THR' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1929 and icode ' ' and resname 'PHE' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1931 and icode ' ' and resname 'LEU' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1932 and icode ' ' and resname 'ASP' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1933 and icode ' ' and resname 'VAL' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1934 and icode ' ' and resname 'ARG' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1935 and icode ' ' and resname 'LEU' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1936 and icode ' ' and resname 'VAL' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1937 and icode ' ' and resname 'PHE' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1938 and icode ' ' and resname 'ASP' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1939 and icode ' ' and resname 'ASN' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1940 and icode ' ' and resname 'CYS' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1941 and icode ' ' and resname 'GLU' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1942 and icode ' ' and resname 'THR' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1943 and icode ' ' and resname 'PHE' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1944 and icode ' ' and resname 'ASN' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1945 and icode ' ' and resname 'GLU' and altid 'A' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1945 and icode ' ' and resname 'GLU' and altid 'B' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1946 and icode ' ' and resname 'ASP' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1947 and icode ' ' and resname 'ASP' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1948 and icode ' ' and resname 'SER' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1949 and icode ' ' and resname 'ASP' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    "chain 'A' and resseq 1950 and icode ' ' and resname 'ILE' and altid ' ' and not (name C or name CA or name N or name O or name CB)",
                    ],
                },
            ],
        ]:

        gs = GenerateLevelSelectionsTask(
            auto_levels = input_levels,
            custom_levels = None,
            overall_selection = None,
            cbeta_in_backbone = True,
            )
        gs.run(hierarchy=hierarchy)

        assert gs.result.level_labels == level_labels
        for i, l in enumerate(level_labels):
            assert gs.result.level_group_selection_strings[i] == sorted(correct_groups[l])

def test_GenerateLevelSelectionsTask_auto_levels_2():

    from pandemic.adp.hierarchy.level_selections import GenerateLevelSelectionsTask

    hierarchy = iotbx.pdb.hierarchy.input(pdb_string=pdb_test_structure_atoms[0]).hierarchy

    for input_levels, level_labels, correct_groups in [
        [
            ['auto_group'],
            ['groups'],
            {
                'groups' : [
                    "chain 'A' and (resid 1900  through 1910 )",
                    "chain 'A' and (resid 1911  through 1950 )",
                    ],
                },
            ],
        ]:

        gs = GenerateLevelSelectionsTask(
            auto_levels = input_levels,
            custom_levels = None,
            overall_selection = None,
            cbeta_in_backbone = True,
            )
        gs.run(hierarchy=hierarchy)

        assert gs.result.level_labels == level_labels
        for i, l in enumerate(level_labels):
            assert gs.result.level_group_selection_strings[i] == sorted(correct_groups[l])

def test_GenerateLevelSelectionsTask_cbeta_inclusion_and_overall_selection():

    from pandemic.adp.hierarchy.level_selections import GenerateLevelSelectionsTask

    hierarchy = iotbx.pdb.hierarchy.input(pdb_string=pdb_test_structure_atoms[0]).hierarchy

    for input_levels, level_labels, correct_groups in [
        [
            ['chain', 'ss', 'residue', 'backbone_sidechain'],
            ['chain', 'sec. struct.', 'residue', 'backbone/sidechain'],
            {
                'chain' : ["chain 'A'"],
                'sec. struct.' : [
                    "chain 'A' and resid 1900  through 1905",
                    "chain 'A' and resid 1906  through 1909",
                    "chain 'A' and resid 1910  through 1920",
                    ],
                'residue' : [
                    "chain 'A' and resseq 1900 and icode ' '",
                    "chain 'A' and resseq 1901 and icode ' '",
                    "chain 'A' and resseq 1902 and icode ' '",
                    "chain 'A' and resseq 1903 and icode ' '",
                    "chain 'A' and resseq 1904 and icode ' '",
                    "chain 'A' and resseq 1905 and icode ' '",
                    "chain 'A' and resseq 1906 and icode ' '",
                    "chain 'A' and resseq 1907 and icode ' '",
                    "chain 'A' and resseq 1908 and icode ' '",
                    "chain 'A' and resseq 1909 and icode ' '",
                    "chain 'A' and resseq 1910 and icode ' '",
                    "chain 'A' and resseq 1911 and icode ' '",
                    "chain 'A' and resseq 1912 and icode ' '",
                    "chain 'A' and resseq 1913 and icode ' '",
                    "chain 'A' and resseq 1914 and icode ' '",
                    "chain 'A' and resseq 1915 and icode ' '",
                    "chain 'A' and resseq 1916 and icode ' '",
                    "chain 'A' and resseq 1917 and icode ' '",
                    "chain 'A' and resseq 1918 and icode ' '",
                    "chain 'A' and resseq 1919 and icode ' '",
                    "chain 'A' and resseq 1920 and icode ' '",
                    ],
                'backbone/sidechain' : [
                    "chain 'A' and resseq 1901 and icode ' ' and resname 'TYR' and altid ' ' and (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1902 and icode ' ' and resname 'LYS' and altid ' ' and (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1903 and icode ' ' and resname 'LYS' and altid ' ' and (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1904 and icode ' ' and resname 'VAL' and altid ' ' and (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1905 and icode ' ' and resname 'ILE' and altid ' ' and (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1906 and icode ' ' and resname 'LYS' and altid ' ' and (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1907 and icode ' ' and resname 'LYS' and altid ' ' and (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1909 and icode ' ' and resname 'MET' and altid ' ' and (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1910 and icode ' ' and resname 'ASP' and altid ' ' and (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1911 and icode ' ' and resname 'PHE' and altid ' ' and (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1912 and icode ' ' and resname 'SER' and altid ' ' and (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1913 and icode ' ' and resname 'THR' and altid ' ' and (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1914 and icode ' ' and resname 'ILE' and altid ' ' and (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1915 and icode ' ' and resname 'ARG' and altid ' ' and (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1916 and icode ' ' and resname 'GLU' and altid ' ' and (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1917 and icode ' ' and resname 'LYS' and altid ' ' and (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1918 and icode ' ' and resname 'LEU' and altid ' ' and (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1919 and icode ' ' and resname 'SER' and altid ' ' and (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1920 and icode ' ' and resname 'SER' and altid ' ' and (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1901 and icode ' ' and resname 'TYR' and altid ' ' and not (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1902 and icode ' ' and resname 'LYS' and altid ' ' and not (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1903 and icode ' ' and resname 'LYS' and altid ' ' and not (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1904 and icode ' ' and resname 'VAL' and altid ' ' and not (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1905 and icode ' ' and resname 'ILE' and altid ' ' and not (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1906 and icode ' ' and resname 'LYS' and altid ' ' and not (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1907 and icode ' ' and resname 'LYS' and altid ' ' and not (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1909 and icode ' ' and resname 'MET' and altid ' ' and not (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1910 and icode ' ' and resname 'ASP' and altid ' ' and not (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1911 and icode ' ' and resname 'PHE' and altid ' ' and not (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1912 and icode ' ' and resname 'SER' and altid ' ' and not (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1913 and icode ' ' and resname 'THR' and altid ' ' and not (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1914 and icode ' ' and resname 'ILE' and altid ' ' and not (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1915 and icode ' ' and resname 'ARG' and altid ' ' and not (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1916 and icode ' ' and resname 'GLU' and altid ' ' and not (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1917 and icode ' ' and resname 'LYS' and altid ' ' and not (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1918 and icode ' ' and resname 'LEU' and altid ' ' and not (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1919 and icode ' ' and resname 'SER' and altid ' ' and not (name C or name CA or name N or name O)",
                    "chain 'A' and resseq 1920 and icode ' ' and resname 'SER' and altid ' ' and not (name C or name CA or name N or name O)",
                    ],
                },
            ],
        ]:

        gs = GenerateLevelSelectionsTask(
            auto_levels = input_levels,
            custom_levels = None,
            overall_selection = 'resseq 1900:1920',
            cbeta_in_backbone = False,
            )
        gs.run(hierarchy=hierarchy)

        assert gs.result.level_labels == level_labels
        for i, l in enumerate(level_labels):
            assert gs.result.level_group_selection_strings[i] == sorted(correct_groups[l])

def test_GenerateLevelSelectionsTask_custom_levels():

    from pandemic.adp.hierarchy.level_selections import GenerateLevelSelectionsTask

    hierarchy = iotbx.pdb.hierarchy.input(pdb_string=pdb_test_structure_atoms[0]).hierarchy

    standard_levels = ['chain', 'residue']
    custom_levels = [
        group_args(
            depth=1,
            insert_before = None,
            insert_after = None,
            label='first_level',
            selection=['resseq 1910', 'resseq 1913'],
            ),
        group_args(
            depth=3,
            insert_before = None,
            insert_after = None,
            label='middle_level',
            selection=['resseq 1911:1920', 'resseq 1921:1930'],
            ),
        group_args(
            depth=5,
            insert_before = None,
            insert_after = None,
            label='last_level',
            selection=['resname ALA', 'resname ILE'],
            ),
    ]

    gs = GenerateLevelSelectionsTask(
        auto_levels = standard_levels,
        custom_levels = custom_levels,
        overall_selection = None,
        cbeta_in_backbone = True,
        )
    gs.run(hierarchy=hierarchy)

    assert gs.result.level_labels == ['first_level', 'chain', 'middle_level', 'residue','last_level']
    assert gs.result.level_group_selection_strings[0] == ['resseq 1910', 'resseq 1913']
    assert gs.result.level_group_selection_strings[1] == ["chain 'A'"]
    assert gs.result.level_group_selection_strings[2] == ['resseq 1911:1920', 'resseq 1921:1930']
    assert len(gs.result.level_group_selection_strings[3]) == 51
    assert gs.result.level_group_selection_strings[4] == ['resname ALA', 'resname ILE']

def test_BuildLevelArrayTask():

    import numpy

    from pandemic.adp.hierarchy.level_array import BuildLevelArrayTask

    hierarchy = iotbx.pdb.hierarchy.input(pdb_string=pdb_test_structure_atoms[0]).hierarchy

    bl = BuildLevelArrayTask(
        overall_selection = None,
        )

    bl.run(
        hierarchy = hierarchy,
        level_group_selection_strings = [
            ['chain A'],
            ['resseq 1911:1920','resseq 1921:1930','resseq 1940:1950','resseq 2000:2010'],
            ['resseq {}'.format(i) for i in range(1900, 1955)],
            ],
        )

    assert bl.result.level_group_array.shape == (3, 428)
    #
    assert (bl.result.level_group_array[0] == 0).all()
    #
    assert set(bl.result.level_group_array[1]) == {-1,0,1,2}
    assert numpy.where(bl.result.level_group_array[1]==-1)[0].tolist() == range(0,90) + range(259,335)
    assert numpy.where(bl.result.level_group_array[1]==+0)[0].tolist() == range(90,171)
    assert numpy.where(bl.result.level_group_array[1]==+1)[0].tolist() == range(171,259)
    assert numpy.where(bl.result.level_group_array[1]==+2)[0].tolist() == range(335,428)
    assert numpy.where(bl.result.level_group_array[1]==+3)[0].tolist() == []
    #
    assert set(bl.result.level_group_array[2]) == set(range(0, 51))
    assert numpy.where(bl.result.level_group_array[2]==-1)[0].tolist() == []
    assert numpy.where(bl.result.level_group_array[2]==+0)[0].tolist() == range(0,4)
    assert numpy.where(bl.result.level_group_array[2]==+1)[0].tolist() == range(4,16)
    assert numpy.where(bl.result.level_group_array[2]==21)[0].tolist() == range(171,175)
    assert numpy.where(bl.result.level_group_array[2]==49)[0].tolist() == range(412,420)
    assert numpy.where(bl.result.level_group_array[2]==50)[0].tolist() == range(420,428)
    assert numpy.where(bl.result.level_group_array[2]==51)[0].tolist() == []

    assert bl.result.overall_atom_mask.shape == (428,)
    assert bl.result.overall_atom_mask.all()

    assert len(bl.result.level_group_selection_strings) == 3
    assert bl.result.level_group_selection_strings[0] == ['chain A']
    assert bl.result.level_group_selection_strings[1] == ['resseq 1911:1920', 'resseq 1921:1930', 'resseq 1940:1950']
    assert bl.result.level_group_selection_strings[2] == ['resseq {}'.format(i) for i in range(1900, 1951)]

    bl = BuildLevelArrayTask(
        overall_selection = 'resseq 1800:1915',
        )

    bl.run(
        hierarchy = hierarchy,
        level_group_selection_strings = [
            ['chain A'],
            ['resseq 1911:1920','resseq 1921:1930','resseq 1940:1950','resseq 2000:2010'],
            ],
        )

    assert bl.result.level_group_array.shape == (2, 133)
    #
    assert (bl.result.level_group_array[0] == 0).all()
    #
    set(bl.result.level_group_array[1]) == {-1,0}
    assert numpy.where(bl.result.level_group_array[1]==-1)[0].tolist() == range(0,90)
    assert numpy.where(bl.result.level_group_array[1]==+0)[0].tolist() == range(90,133)

    assert bl.result.overall_atom_mask.shape == (428,)
    assert numpy.where(bl.result.overall_atom_mask)[0].tolist() == range(0,133)

    assert len(bl.result.level_group_selection_strings) == 2
    assert bl.result.level_group_selection_strings[0] == ['chain A']
    assert bl.result.level_group_selection_strings[1] == ['resseq 1911:1920']

    bl = BuildLevelArrayTask(
        overall_selection = 'resseq 1800:1915',
        )

    with raises(Exception) as e:
        bl.run(
            hierarchy = hierarchy,
            level_group_selection_strings = [
                ['chain A'],
                ['resseq 1920:1930'],
                ],
            )
    assert str(e.value) == "Levels have been created that contain no atoms!\nInvalid Levels: [2]"

def test_BuildLevelArrayAsTreeTask():

    import numpy

    from pandemic.adp.hierarchy.level_array_tree import BuildLevelArrayAsTreeTask

    level_group_array = -1 * numpy.ones((3, 20))

    level_group_array[0, 00:15] = 0

    level_group_array[1, 00:05] = 0
    level_group_array[1, 10:15] = 1

    level_group_array[2, 00:02] = 0
    level_group_array[2, 02:05] = 1
    level_group_array[2, 05:07] = 2
    level_group_array[2, 07:10] = 3
    level_group_array[2, 10:13] = 4
    level_group_array[2, 13:15] = 5
    level_group_array[2, 15:18] = 6

    bl = BuildLevelArrayAsTreeTask()

    bl.run(level_group_array = level_group_array)

    assert bl.result.tree.links == {
        0: {0: {1: [0, 1], 2: [2, 3]}},
        1: {0: {2: [0, 1]}, 1: {2: [4, 5]}},
        2: {0: {}, 1: {}, 2: {}, 3: {}, 4: {}, 5: {}, 6: {}},
    }

def test_CreateHierarchicalModelTask():
    """Integration test"""

    import numpy

    from pandemic.adp.hierarchy import CreateHierarchicalModelTask

    hierarchy = iotbx.pdb.hierarchy.input(pdb_string=pdb_test_structure_atoms[0]).hierarchy

    hm = CreateHierarchicalModelTask(
        auto_levels = ['chain', 'ss'],
        custom_levels = None,
        overall_selection = 'resseq 1900:1920',
        cbeta_in_backbone = True,
        )

    hm.run(hierarchy)

    assert hm.result.level_labels == ['chain', 'sec. struct.']
    #
    assert hm.result.level_group_array.shape == (2, 171)
    #
    assert (hm.result.level_group_array[0] == 0).all()
    #
    assert set(hm.result.level_group_array[1]) == {0,1,2}
    assert numpy.where(hm.result.level_group_array[1]==0)[0].tolist() == range(0,49)
    assert numpy.where(hm.result.level_group_array[1]==1)[0].tolist() == range(49,82)
    assert numpy.where(hm.result.level_group_array[1]==2)[0].tolist() == range(82,171)
    #
    assert len(hm.result.level_group_selection_strings) == 2
    assert hm.result.level_group_selection_strings[0] == ["chain 'A'"]
    assert hm.result.level_group_selection_strings[1] == ["chain 'A' and resid 1900  through 1905", "chain 'A' and resid 1906  through 1909", "chain 'A' and resid 1910  through 1920"]
    #
    assert hm.result.level_group_tree.links == {0: {0: {1: [0, 1, 2]}}, 1: {0: {}, 1: {}, 2: {}}}
    #
    assert hm.result.overall_atom_mask.shape == (428,)
    assert numpy.where(hm.result.overall_atom_mask)[0].tolist() == range(0,171)

def test_translate_phenix_selections_to_pymol_selections_simple():

    from pandemic.adp.hierarchy.summary import translate_phenix_selections_to_pymol_selections_simple

    for phenix_str, pymol_str in [
            (
                ['chain B'],
                ['chain B'],
                ),
            (
                ['chain C and resid 10 through 21'],
                ['chain C and resi 10:21'],
                ),
            (
                ['chain B and resid 121B through 130'],
                ['chain B and resi 121B:130'],
                ),
            (
                ['chain A and resseq 15:20'],
                [None],
                ),
            (
                ['chain A', 'chains B'],
                ['chain A', None],
                ),
            ]:

        t = translate_phenix_selections_to_pymol_selections_simple(phenix_str)
        assert t == pymol_str

