import shutil
import pathlib as pl


class MultiStateResources(object):

    def __init__(self, dir_path):

        self.dir_path = pl.Path(str(dir_path))

    def get(self, rel_filepath, out_path=None):

        in_path = (
            self.dir_path / rel_filepath
            )

        assert in_path.exists()

        if out_path is not None:
            shutil.copy(
                str(in_path),
                str(out_path),
                )

            return out_path

        return in_path

    def read(self, path):

        with open(str(path), 'r') as fh: 
            s = fh.read()

        return s.strip()

    def assert_same(self, path1, path2):

        string1 = self.read(path1)
        string2 = self.read(path2)

        assert string1 == string2, "{} != {}".format(
            str(path1), str(path2),
            )

    ###

    def input_apo(self, out_path=None):
        return self.get("states/model-apo.pdb", out_path)

    def input_glc(self, out_path=None):
        return self.get("states/model-GLC.pdb", out_path)

    def input_glo(self, out_path=None):
        return self.get("states/model-GLO.pdb", out_path)

    def input_z9n(self, out_path=None):
        return self.get("states/model-Z9N.pdb", out_path)

    ###

    def merged_apo_glc(self, out_path=None):
        return self.get("merged/apo_glc.pdb", out_path)

    def merged_apo_glc_glo(self, out_path=None):
        return self.get("merged/apo_glc_glo.pdb", out_path)

    def merged_apo_glc_glo_z9n(self, out_path=None):
        return self.get("merged/apo_glc_glo_z9n.pdb", out_path)

    ###

    def restraints_apo_glc_phenix(self, out_path=None):
        return self.get("restraints/apo_glc.phenix.params", out_path)

    def restraints_apo_glc_refmac(self, out_path=None):
        return self.get("restraints/apo_glc.refmac.params", out_path)

    def restraints_apo_glc_glo_phenix(self, out_path=None):
        return self.get("restraints/apo_glc_glo.phenix.params", out_path)

    def restraints_apo_glc_glo_refmac(self, out_path=None):
        return self.get("restraints/apo_glc_glo.refmac.params", out_path)

    def restraints_apo_glc_glo_z9n_phenix(self, out_path=None):
        return self.get("restraints/apo_glc_glo_z9n.phenix.params", out_path)

    def restraints_apo_glc_glo_z9n_refmac(self, out_path=None):
        return self.get("restraints/apo_glc_glo_z9n.refmac.params", out_path)

    #

    ###

    def split_apo_glc_to_apo(self, out_path=None):
        return self.get("split/apo_glc-split-ABC.pdb", out_path)

    def split_apo_glc_to_glc(self, out_path=None):
        return self.get("split/apo_glc-split-DEF.pdb", out_path)

    def split_apo_glc_to_apo_glc(self, out_path=None):
        return self.get("split/apo_glc-split-ABCDEF.pdb", out_path)

    def split_apo_glc_to_A(self, out_path=None):
        return self.get("split/apo_glc-split-A.pdb", out_path)

    def split_apo_glc_to_B(self, out_path=None):
        return self.get("split/apo_glc-split-B.pdb", out_path)

    def split_apo_glc_to_C(self, out_path=None):
        return self.get("split/apo_glc-split-C.pdb", out_path)

    def split_apo_glc_to_D(self, out_path=None):
        return self.get("split/apo_glc-split-D.pdb", out_path)

    def split_apo_glc_to_E(self, out_path=None):
        return self.get("split/apo_glc-split-E.pdb", out_path)

    def split_apo_glc_to_F(self, out_path=None):
        return self.get("split/apo_glc-split-F.pdb", out_path)

    #

    def split_apo_glc_to_apo_no_prune(self, out_path=None):
        return self.get("split/apo_glc_no_prune-split-ground-state-ABC.pdb", out_path)

    def split_apo_glc_to_glc_no_prune(self, out_path=None):
        return self.get("split/apo_glc_no_prune-split-bound-state-DEF.pdb", out_path)

    def split_apo_glc_to_apo_no_reset_altlocs(self, out_path=None):
        return self.get("split/apo_glc_no_reset_altlocs-split-ground-state-ABC.pdb", out_path)

    def split_apo_glc_to_glc_no_reset_altlocs(self, out_path=None):
        return self.get("split/apo_glc_no_reset_altlocs-split-bound-state-DEF.pdb", out_path)

    #

    def split_apo_glc_glo_to_apo(self, out_path=None):
        return self.get("split/apo_glc_glo-split-ABC.pdb", out_path)

    def split_apo_glc_glo_to_glc(self, out_path=None):
        return self.get("split/apo_glc_glo-split-DEF.pdb", out_path)

    def split_apo_glc_glo_to_glo(self, out_path=None):
        return self.get("split/apo_glc_glo-split-GHI.pdb", out_path)

    def split_apo_glc_glo_to_apo_glc(self, out_path=None):
        return self.get("split/apo_glc_glo-split-ABCDEF.pdb", out_path)

    def split_apo_glc_glo_to_glc_glo(self, out_path=None):
        return self.get("split/apo_glc_glo-split-DEFGHI.pdb", out_path)


