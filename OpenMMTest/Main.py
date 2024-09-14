# Import the needed libraries
import sys
import os
import random

from openmm import *
from openmm import app
from openmm import unit
from pdbfixer import PDBFixer

def main():
    # Load a PDB file from the DataPreparation/PDBData folder
    folder_path = 'DataPreparation/PDBData'

    # Randomly select a file from the list
    random_file = random.choice(os.listdir(folder_path))

    # Construct the full file path
    pdb_file = os.path.join(folder_path, random_file)
    
    # Fix the PDB file
    fixer = PDBFixer(filename=pdb_file)
    
    fixer.findMissingResidues()
    
    # TODO: Some force fields do not recognize certain residues. Will have to find a workaround for this
    # Replace residues with the most common residue, in this case ALA
    fixer.findNonstandardResidues()
    fixer.nonstandardResidues = [(residue, 'ALA') for residue, replacement in fixer.nonstandardResidues]
    fixer.replaceNonstandardResidues()
    
    fixer.removeHeterogens()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens()
    
    # Create a folder to store the fixed PDB files
    if not os.path.exists('OpenMMTest/FixedPDBData'):
        os.makedirs('OpenMMTest/FixedPDBData')
    
    # Save the fixed PDB file
    fixed_pdb_file = 'OpenMMTest/FixedPDBData/fixed_' + random_file
    app.PDBFile.writeFile(fixer.topology, fixer.positions, open(fixed_pdb_file, 'w'))

    # Load the PDB file
    pdb = app.PDBFile(fixed_pdb_file) 
    
    # NOTE: For now go with a set force field. Ask the user what type of force field they want to use
    # force_field = input("Which force field would you like to use? (amber14-all.xml, charmm36.xml, or oplsaa.xml): ")
    force_field = 'amber14-all.xml'

    # NOTE: For now go with a set water model. Ask the user what type of water model they want to use
    # water_model = input("Which water model would you like to use? (tip3pfb.xml, tip4pew.xml, or spce.xml): ")
    water_model = 'amber14/tip3pfb.xml'
    
    # Create a force field, first file is the force field, second file is the water model
    forcefield = app.ForceField(force_field, water_model)

    # Create the system
    system = forcefield.createSystem(topology=pdb.topology, 
                                     nonbondedMethod=app.PME, 
                                     nonbondedCutoff=unit.Quantity(value=1.0, unit=unit.nano*unit.meters))

    # Set up the integrator
    integrator = LangevinIntegrator(300*unit.kelvin, 1/(unit.pico*unit.seconds), 0.002*(unit.pico*unit.seconds))

    # Create the simulation
    simulation = app.Simulation(pdb.topology, system, integrator)

    # Set the initial positions
    simulation.context.setPositions(pdb.positions)

    # Minimize energy
    print("Minimizing energy...")
    simulation.minimizeEnergy()

    # Set up the reporter to save the trajectory
    simulation.reporters.append(app.PDBReporter('OpenMMTest/output.pdb', 100))
    simulation.reporters.append(app.StateDataReporter(sys.stdout, 100, step=True, potentialEnergy=True, temperature=True))

    # Run the simulation
    print("Running simulation...")
    simulation.step(1000)  # Adjust the number of steps as needed

if __name__ == "__main__":
    main()