import sys
from rdkit import Chem

# Use the actual AlphaFold PDB file as the template
template_file = 'AF-Q43127-F1.pdb'
mol = Chem.MolFromPDBFile(template_file, removeHs=False)

if mol is None:
    print(f"Error: Could not find {template_file}. Ensure it's in the current directory.")
    sys.exit(1)

# In the wild-type, residue 125 is naturally Aspartic Acid (ASP)
# We ensure the metadata reflects this for the baseline
for atom in mol.GetAtoms():
    res_info = atom.GetPDBResidueInfo()
    if res_info is not None and res_info.GetResidueNumber() == 125:
        res_info.SetResidueName("ASP")

with open('gs2_wildtype.pdb', 'w') as f:
    f.write(Chem.MolToPDBBlock(mol))
    f.flush()

print("Successfully created gs2_wildtype.pdb using AlphaFold template.")