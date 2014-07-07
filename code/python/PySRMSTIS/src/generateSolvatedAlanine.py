
#!/usr/local/bin/env python

"""
Part of a tutorial program illustrating how to set up a simple MD simulation with OpenMM

DESCRIPTION
solvating a pdb structure
@author Antonia Mey <antonia.mey@fu-berlin.de>

"""



#imports
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import time

numPlatforms = Platform.getNumPlatforms()
print("There are", numPlatforms, "Platforms available:")
print()

platform = openmm.Platform.getPlatformByName('CPU') # platform to use <- change accordingly to Cuda or CPU
print ('Platform used is ' + platform.getName())

start_time = time.time()

#loading the unsolvated pdb file
print('Loading...')
pdb = PDBFile('data/input.pdb')
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
modeller = Modeller(pdb.topology, pdb.positions)

#choosing the solvent
print('Adding solvent...')
modeller.addSolvent(forcefield, model='tip3p', padding=1*nanometer)
print('Minimizing...')
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME)
integrator = VerletIntegrator(0.001*picoseconds)
simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy(maxIterations=200)
#saving the solvated pdb file
print('Saving...')
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('data/Alanine_solvated.pdb', 'w'))

end_time = time.time()
print "total time "+str(end_time-start_time)	
