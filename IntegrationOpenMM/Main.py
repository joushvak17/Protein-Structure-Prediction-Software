from openmm.app import *
from openmm.unit import nano, meters, kelvin, pico
from openmm import *
from sys import stdout

pdb = PDBFile("PDBTest.pdb")
forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
                                 nonbondedCutoff=1.0*nano*meters,
                                 constraints=HBonds)
integrator = LangevinMiddleIntegrator(300.0*kelvin, 1.0/pico,
                                      0.004*pico)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter("output.pdb",
                                        1000))
simulation.reporters.append(StateDataReporter(stdout, 1000,
                                              step=True,
                                              potentialEnergy=True,
                                              temperature=True))
simulation.step(10000)
