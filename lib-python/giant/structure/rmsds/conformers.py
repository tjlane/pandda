import giant.logs as lg
logger = lg.getLogger(__name__)

import numpy as np

from .atoms import (
    CalculatePairedAtomSetsRMSD,
    )

class CalculatePairedConformerSetsRMSDs(object):

    def __init__(self):
        
        self.calculate_paired_atom_sets_rmsd = (
            CalculatePairedAtomSetsRMSD()
            )

    def __call__(self,
        conformers_1,
        conformers_2 = None,
        ):

        if conformers_2 is None:
            conformers_2 = conformers_1

        rmsds = np.empty(
            (len(conformers_1), len(conformers_2)),
            dtype = float,
            )

        rmsds[:] = None

        for i, c1 in enumerate(conformers_1):

            for j, c2 in enumerate(conformers_2):
                
                rmsds[i,j] = self.calculate_paired_atom_sets_rmsd(
                    atoms_1 = c1.atoms(),
                    atoms_2 = c2.atoms(),
                    )

        ret_list = []

        while not np.isnan(rmsds).all():

            min_val = np.nanmin(rmsds)

            i1, i2 = list(zip(*
                np.where(rmsds==min_val)
                ))[0]

            # Clear these values so that a conformer cannot be used more than once
            rmsds[i1, :] = None
            rmsds[:, i2] = None

            # Return conformers and the rmsd
            ret_list.append(
                (conformers_1[i1], conformers_2[i2], min_val)
                )

            # logger(rmsds)
            # logger(ret_list[-1])

        return ret_list
