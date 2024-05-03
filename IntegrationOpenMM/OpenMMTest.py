# Import the needed libraries
import pandas as pd

from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

# Input the test PDB file to a PDBFile object
pdb = PDBFile("PDBFix.pdb")

# Read the dataset file and turn it into a dictionary
pdb_data = pd.read_csv("PDBOpenMMTest.csv").to_dict(orient="records")

# Loop through the dictionary
for entry in pdb_data:
    # Check to see if the experimental method matches
    if entry["Experimental"] == "x-ray diffraction":
        # Input the force fields that will be used, in this case AMBER14 for protein
        forcefield = ForceField("", "amber14/tip3pfb.xml")
    elif entry["Experimental"] == "neutron diffraction":
        # Input the alternative AMBER 14 force field which is provided by OpenMM
        forcefield = ForceField("amber14/protein.ff15ipq.xml", "amber14/tip4pfb.xml")

    # Decide on the method to use for nonbonded methods based on resolution
    if entry["Resolution"] <= 1.89:
        nonbonded_method = CutoffPeriodic
    elif entry["Resolution"] >= 2.60:
        nonbonded_method = LJPME
    else:
        nonbonded_method = PME

# Combine the force field with the molecular topology from the PDB file
system = forcefield.createSystem(pdb.topology,
                                 nonbondedMethod=nonbonded_method,
                                 nonbondedCutoff=Quantity(value=1.2, unit=nano*meter),
                                 constraints=HBonds)

# Create an integrator to use for the equations of motion
integrator = LangevinMiddleIntegrator(300*kelvin,
                                      1/Quantity(value=1.0, unit=pico*second),
                                      Quantity(value=0.002, unit=pico*second))

# Combine everything so a Simulation object is created
simulation = Simulation(pdb.topology, system, integrator)

# Specify the initial atom positions for the simulation
simulation.context.setPositions(pdb.positions)

# Perform local energy minimization
simulation.minimizeEnergy()

# Generate the output during from the simulation
simulation.reporters.append(PDBReporter("output.pdb",
                                        1000))

# Get regular status updates from the simulation
simulation.reporters.append(StateDataReporter(stdout,
                                              1000,
                                              step=True,
                                              potentialEnergy=True,
                                              temperature=True))

# Run the simulation
simulation.step(10000)
