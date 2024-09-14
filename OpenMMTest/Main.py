# Import the needed libraries
import sys
import os
import random

from openmm import *
from openmm import app
from openmm import unit

def main():
    # Load a PDB file from the DataPreparation/PDBData folder
    folder_path = 'DataPreparation/PDBData'

    # Randomly select a file from the list
    random_file = random.choice(os.listdir(folder_path))

    # Construct the full file path
    pdb_file = os.path.join(folder_path, random_file)

    # Load the PDB file
    pdb = app.PDBFile(pdb_file) 
    
    # NOTE: For now go with a set force field. Ask the user what type of force field they want to use
    # force_field = input("Which force field would you like to use? (amber14-all.xml, charmm36.xml, or oplsaa.xml): ")
    force_field = 'amber14-all.xml'

    # NOTE: For no go with a set water model. Ask the user what type of water model they want to use
    # water_model = input("Which water model would you like to use? (tip3pfb.xml, tip4pew.xml, or spce.xml): ")
    water_model = 'amber14/tip3pfb.xml'
    
    # Create a force field, first file is the force field, second file is the water model
    forcefield = app.ForceField(force_field, water_model)

    # TODO: Will need to fix the PDB files first
    # Create the system
    system = forcefield.createSystem(topology=pdb.topology, nonbondedMethod=app.PME, nonbondedCutoff=unit.Quantity(value=1.0, unit=unit.nano*unit.meters))

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