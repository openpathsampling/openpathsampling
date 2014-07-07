'''
Simple OpenMM Example

Loading PDB, Solvating, Minimizing and short production run

@author: Jan-Hendrik Prinz
@author: Antonia Mey
'''


#=============================================================================================
# Import OpenMM
#=============================================================================================

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

#=============================================================================================
# Misc Imports
#=============================================================================================

import time
from sys import stdout


start_time = time.time()

#=============================================================================================
# Available OpenMM Platforms
#=============================================================================================

platform = openmm.Platform.getPlatformByName('CPU') # platform to use <- change accordingly to 'Cuda' or 'CPU' or 'Reference' (VERY slow)
print ('Platform used is : ' + platform.getName())


#=============================================================================================
# Loading the unsolvated pdb file
#=============================================================================================

print('Loading...')
pdb = PDBFile('data/input.pdb')
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')


#=============================================================================================
# Adding solvent
#=============================================================================================

print('Adding solvent...')
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(forcefield, model='tip3p', padding=1*nanometer)

#=============================================================================================
# Create System
#=============================================================================================

system1 = forcefield.createSystem(modeller.topology, nonbondedMethod=PME)


#=============================================================================================
# Minimizing
#=============================================================================================

print('Minimizing using VerletIntegrator ...')
integrator1 = VerletIntegrator(0.001*picoseconds)
simulation1 = Simulation(modeller.topology, system1, integrator1, platform)
simulation1.context.setPositions(modeller.positions)
simulation1.minimizeEnergy(maxIterations=200)


#=============================================================================================
# Simulate
#=============================================================================================

print('Simulate using LangevinIntegrator ...')
system2 = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
integrator2 = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)


#=============================================================================================
# Creating simulation context
#=============================================================================================

simulation2 = Simulation(modeller.topology, system2, integrator2) 
simulation2.context.setPositions(simulation1.context.getState(getPositions=True).getPositions())


#=============================================================================================
# Minimizing System 
#=============================================================================================

simulation2.minimizeEnergy(maxIterations=25)


#=============================================================================================
# Adding Reporters 
#=============================================================================================

simulation2.reporters.append(PDBReporter('output_exercise1.pdb', 5)) 
simulation2.reporters.append(StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True))


#=============================================================================================
# Running simulation
#=============================================================================================

simulation2.step(1000)

end_time = time.time()