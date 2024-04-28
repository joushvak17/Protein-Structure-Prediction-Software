from pdbfixer import PDBFixer
from openmm.app import PDBFile

# Inpur the PDB file into the PDBFixer object
fixer = PDBFixer(filename="PDBTest.pdb")

# Find the missing residues
fixer.findMissingResidues()

# Find and replace the nonstandard residues
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()

# Remove the heterogens, but keep the water molecules
fixer.removeHeterogens(True)

# Find and add missing atoms
fixer.findMissingAtoms()
fixer.addMissingAtoms()

# Add missing hydrogens and pass 7.0 for the pH
fixer.addMissingHydrogens(7.0)

# Write the fixed PDB file
PDBFile.writeFile(fixer.topology, fixer.positions, open("PDBFix.pdb", "w"))
