# Import the needed libraries
import sys
import os
import random

from openmm import *
from openmm import app

def main():
    # Load a PDB file from the DataPreparation/PDBData folder
    folder_path = 'DataPreparation/PDBData'

    # Randomly select a file from the list
    random_file = random.choice(os.listdir(folder_path))

    # Construct the full file path
    pdb_file = os.path.join(folder_path, random_file)

    # Load the PDB file
    pdb = app.PDBFile(pdb_file) 
    
    # Create a force field
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    # Create the system
    system = forcefield.createSystem(pdb.topology, nonbondedCutoff=1.0*nanometers)

    # Set up the integrator
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

    # Create the simulation
    simulation = app.Simulation(pdb.topology, system, integrator)

    # Set the initial positions
    simulation.context.setPositions(pdb.positions)

    # Minimize energy
    print("Minimizing energy...")
    simulation.minimizeEnergy()

    # Set up the reporter to save the trajectory
    simulation.reporters.append(app.DCDReporter('output.dcd', 1000))
    simulation.reporters.append(app.StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True))

    # Run the simulation
    print("Running simulation...")
    simulation.step(10000)  # Adjust the number of steps as needed

if __name__ == "__main__":
    main()