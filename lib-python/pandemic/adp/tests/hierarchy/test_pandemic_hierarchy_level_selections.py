from pytest import approx, raises, mark

from libtbx import group_args
import iotbx.pdb

from giant.paths import is_available

from pandemic.resources.structure_test_snippets import pdb_test_structure_atoms

only_if_phenix_is_available = mark.skipif(
    not is_available('phenix.find_tls_groups'),
    reason = "phenix.find_tls_groups is not available",
)

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

@only_if_phenix_is_available
def test_GenerateLevelSelectionsTask_auto_levels_2():

    from pandemic.adp.hierarchy.level_selections import GenerateLevelSelectionsTask

    hierarchy = iotbx.pdb.hierarchy.input(pdb_string=pdb_test_structure_atoms[0]).hierarchy

    for input_levels, level_labels, correct_groups in [
        [
            ['phenix_find_tls_groups'],
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

