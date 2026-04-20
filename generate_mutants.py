import os
from Bio.PDB import PDBParser, PDBIO

# Top 5 Candidates from your Research Archive
MUTATIONS = [
    (112, 'ASN', 'TYR', 'N112Y'),
    (117, 'SER', 'CYS', 'S117C'),
    (170, 'GLU', 'ALA', 'E170A'),
    (125, 'ASP', 'ASN', 'D125N'),
    (111, 'TRP', 'GLY', 'W111G')
]

def create_mutant(input_pdb, res_pos, mut_name):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("WT", input_pdb)
    
    # Simple residue renaming for visualization
    # Note: This changes the label, not the atom coordinates
    for residue in structure.get_residues():
        if residue.get_id()[1] == res_pos:
            residue.resname = mut_name[res_pos > 99 and -1 or -3:] # Keep 3-letter code
            
    io = PDBIO()
    io.set_structure(structure)
    io.save(f"gs2_model_{mut_name}.pdb")
    print(f"✅ Created: gs2_model_{mut_name}.pdb")

if __name__ == "__main__":
    template = "AF-Q43127-F1.pdb"
    if os.path.exists(template):
        for pos, wt, mut_aa, label in MUTATIONS:
            create_mutant(template, pos, label)
    else:
        print(f"❌ Error: {template} not found. Run the curl command first.")