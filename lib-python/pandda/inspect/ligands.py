import os, shutil, string
import pathlib as pl

from giant.refinement.restraints.acedrg import (
    generate_restraints,
    )


class MakeNewLigand(object):

    def __init__(self,
        output_directory,
        ):

        self.output_directory = pl.Path(
            output_directory
            )

        if not self.output_directory.exists():
            self.output_directory.mkdir(parents=True)

        self.tmp_directory = (
            self.output_directory / 'make_ligand_tmp'
            )

        self.disallowed = set(string.punctuation)
        self.disallowed.add(' ')

    def __call__(self,
        ligand_id,
        ligand_name,
        ligand_smiles,
        ):

        self.validate(
            ligand_id = ligand_id,
            ligand_name = ligand_name,
            ligand_smiles = ligand_smiles,
            )

        real_out_pdb = (
            self.output_directory / (ligand_name + '.pdb')
            )

        real_out_cif = real_out_pdb.with_suffix('.cif')

        if real_out_pdb.exists() or real_out_cif.exists():
            raise ValueError(
                'Output ligand files already exist -- please choose a different id/name'
                )

        # Remove temporary directory if it exists
        self.cleanup()

        # This should really always evaluate to True...
        if not self.tmp_directory.exists():
            self.tmp_directory.mkdir(parents=True)

        try:

            tmp_pdb, tmp_cif = generate_restraints(
                smiles = ligand_smiles,
                name = ligand_id, 
                prefix = str(
                    self.tmp_directory / ligand_name
                    ),
                )

        except Exception as e: 
            raise Exception(
                'Error during ligand generation. See terminal output for acedrg error message.\n{}'.format(str(e))
                )

        self.move_output_files(
            ligand_pdb_tmp = tmp_pdb,
            ligand_pdb = real_out_pdb,
            ligand_cif_tmp = tmp_cif,
            ligand_cif = real_out_cif,
            )

        # Only do this upon successful completion
        self.cleanup()

        return {
            'pdb' : str(real_out_pdb),
            'cif' : str(real_out_cif),
            'id' : ligand_id,
            'name' : ligand_name,
            'smiles' : ligand_smiles,
        }

    def validate(self, 
        ligand_id,
        ligand_name,
        ligand_smiles,
        ):

        allowed_name = ['-','_']
        disallowed_name = self.disallowed.difference(allowed_name)

        allowed_id = ['_']
        disallowed_id = self.disallowed.difference(allowed_id)

        if disallowed_name.intersection(ligand_name):

            raise ValueError(
                'Ligand name cannot contain spaces or punctuation except for {}'.format(
                    ' or '.join(allowed_path),
                    )
                )

        if disallowed_id.intersection(ligand_id):
            
            raise ValueError(
                'Ligand ID cannot contain spaces or punctuation except for {}'.format(
                    ' or '.join(allowed_id),
                    )
                )

        if len(ligand_id) == 0:

            raise ValueError(
                'No ligand id provided'
                )

        if len(ligand_id) != 3:

            raise ValueError(
                'Ligand ID must be three characters'
                )

        if len(ligand_name) == 0:

            raise ValueError(
                'No ligand name provided'
                )
            
        if len(ligand_smiles) == 0:

            raise ValueError(
                'No ligand smiles provided'
                )
            
        return

    def move_output_files(self,
        ligand_pdb_tmp,
        ligand_pdb,
        ligand_cif_tmp,
        ligand_cif,
        ):

        ligand_pdb_tmp = pl.Path(ligand_pdb_tmp)
        ligand_pdb_tmp.rename(ligand_pdb)

        ligand_cif_tmp = pl.Path(ligand_cif_tmp)
        ligand_cif_tmp.rename(ligand_cif)

        assert ligand_pdb.exists()
        assert ligand_cif.exists()

    def cleanup(self):

        if self.tmp_directory.exists():
            shutil.rmtree(str(self.tmp_directory))


