# Import the needed libraries
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

# Input the test PDB file to a PDBFile object
pdb = PDBFile("PDBFix.pdb")

# Input the force fields that will be used, in this case AMBER14
forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

# Combine the force field with the molecular topology from the PDB file
system = forcefield.createSystem(pdb.topology,
                                 nonbondedMethod=PME,
                                 nonbondedCutoff=None,
                                 constraints=HBonds)

# Create an integrator to use for the equations of motion
integrator = LangevinMiddleIntegrator(300 * kelvin,
                                      1 / pico,
                                      0.002 * pico)

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
