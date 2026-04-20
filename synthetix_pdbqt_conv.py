import os
from meeko import MoleculePreparation
from rdkit import Chem

def prepare_ligand(pdb_file, output_pdbqt):
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
    prepper = MoleculePreparation()
    prepper.prepare(mol)
    prepper.write_pdbqt_file(output_pdbqt)
    print(f"Prepared ligand: {output_pdbqt}")

# Prepare the substrates
prepare_ligand("ATP_ligand.pdb", "ATP.pdbqt")
prepare_ligand("Glutamate_ligand.pdb", "Glutamate.pdbqt")

# For the receptor, we use a simple command-line conversion 
# or a specialized tool since it is a static 'rigid' body.
print("Ligands converted. Ready for AutoDock Vina.")