from pdbfixer import PDBFixer
from openmm.app import PDBFile

# Load the validated D125N mutant created on your MacBook Pro M3
# This is a key lead for the Synthetix/Veron project
fixer = PDBFixer(filename='gs2_model_D125N.pdb')

# Identify and repair structural gaps
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()

# Add hydrogens for the alkaline chloroplastic environment
# Note: 'pH' must be uppercase to avoid the TypeError
fixer.addMissingHydrogens(pH=8.0) 

# Save the finalized, physics-ready receptor
with open('D125N_ready.pdb', 'w') as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)

print("D125N model is now 'ready' for docking simulation.")