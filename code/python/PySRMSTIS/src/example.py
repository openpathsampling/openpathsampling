'''
Created on 05.06.2014

@author: Jan-Hendrik Prinz
'''

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout


# Reading structure and force field files
pdb = PDBFile('data/input.pdb')

forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
# Creating System
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
# Creating simulation context
simulation = Simulation(pdb.topology, system, integrator) 

simulation.context.setPositions(pdb.positions)
# Minimizing System 
simulation.minimizeEnergy(maxIterations=25)
# Adding Reporters 
simulation.reporters.append(PDBReporter('output_exercise1.pdb', 5)) 
simulation.reporters.append(StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True))
# Running simulation
simulation.step(1000)
