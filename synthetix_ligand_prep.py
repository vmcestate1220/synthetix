from rdkit import Chem
from rdkit.Chem import AllChem

# Define substrates as SMILES strings
substrates = {
    "ATP": "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
    "Glutamate": "C(CC(=O)O)[C@@H](C(=O)O)N"
}

for name, smiles in substrates.items():
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    Chem.MolToPDBFile(mol, f"{name}_ligand.pdb")
    print(f"Generated 3D structure for: {name}")