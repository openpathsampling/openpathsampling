#!/usr/local/bin/env python

"""
Tutorial program illustrating how to set up a simple MD simulation with OpenMM

DESCRIPTION
Crude REMD simulation with method mixReplicas(iteration)

@author Antonia Mey <antonia.mey@fu-berlin.de>



"""

#=============================================================================================
# PARAMETERS
#=============================================================================================

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import numpy as np
import time
import math

#=============================================================================================
# PARAMETERS
#=============================================================================================
kB = BOLTZMANN_CONSTANT_kB * AVOGADRO_CONSTANT_NA
temperatures = [300.0, 306.7, 313.6, 320.6, 327.9, 335.2, 342.8, 350.5, 358.3, 366.4, 374.6, 383.0, 391.6, 400.4, 409.4, 418.6, 428.0, 437.6, 447.4, 457.4, 467.8, 478.3, 489.0, 500.0]
collision_rate = 1.0 / picoseconds # collision rate for Langevin dynamics
timestep = 2.0 * femtoseconds # timestep for Langevin integrator
nequib_steps = 5 #number of nvt equilibration steps with position constraints on Alanine
#platform = openmm.Platform.getPlatformByName("Cuda") # platform to use

max_iterations = 5 # number of iterations of nsteps_per_interation of dynamics to run
nsteps_per_iteration = 500 # number of steps of dynamics per iteration
Alanine_atoms = 22

#=============================================================================================
# Functions
#=============================================================================================

def mixReplicas(iteration):
	"""  
	neighbouring replicas are being exchanged
	"""
	attempts = 0.0
	accepted = 0.0
        start_time = time.time()
	print "attempting exchanges now.."
	for npair in range(len(simulations)-1):
		attempts=attempts+1
			
		n1 = npair
		n2 = npair+1
		if iteration%2!=0:
			n1 = len(temperatures)-npair-1
			n2 = len(temperatures)-npair-2
		#print "n2 is " + str(n2)
		#print "temperature of n2 is "+ str(integrators[n2].getTemperature())
		#get the potential energy
		ePot1 = simulations[n1].context.getState(getEnergy=True).getPotentialEnergy()
		ePot2 = simulations[n2].context.getState(getEnergy=True).getPotentialEnergy()
		delta = ((1.0 / (integrators[n1].getTemperature()*kB)) - (1.0 / (integrators[n2].getTemperature()*kB)))*(ePot1 - ePot2)
		if (np.random.rand() <math.exp(1.0*delta)):
			accepted = accepted+1
			temp = integrators[n1].getTemperature();
                	integrators[n1].setTemperature(integrators[n2].getTemperature());
                	integrators[n2].setTemperature(temp);
	print "Percentage of accepted exchanges "+str((accepted/attempts)*100.0)
	#reassign all the velocities
	print "reassinging velocities"
	for replica in range(len(simulations)):
  		nparticles = system.getNumParticles()

		temperature = integrators[replica].getTemperature()
        	kT = kB * temperature # thermal energy    
        	sqrt_kT_over_m = Quantity(numpy.zeros([nparticles,3], numpy.float64), nanometers / picosecond)
        	for atom_index in range(nparticles):
            		mass = system.getParticleMass(atom_index) # atomic mass
            		sqrt_kT_over_m[atom_index,:] = sqrt(kT / mass) # standard deviation of velocity distribution for each coordinate for this atom
		velocities = sqrt_kT_over_m * numpy.random.standard_normal(size=(nparticles,3))
		#print velocities
		simulations[replica].context.setVelocities(velocities)  




stime = time.time()
#equilibration at different temperatures
pdb = PDBFile('Alanine_solvated.pdb')
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
integrators = []
for i in range (len(temperatures)):
	integrators.append(LangevinIntegrator(temperatures[i], collision_rate, timestep))

simulations = []
for i in range(len(temperatures)):
	simulations.append(Simulation(pdb.topology, system, integrators[i]))
	simulations[i].context.setPositions(pdb.positions)

#Dirty Position Constraints for NVT
Alanine_masses = np.zeros(Alanine_atoms, numpy.double)
Alanine_masses = Quantity(Alanine_masses, dalton)
for i in range(Alanine_atoms):
	Alanine_masses[i]= system.getParticleMass(i)
	system.setParticleMass(i, 0.0)
for i in range(len(simulations)):
	print "We are at temperature " + str(temperatures[i])
	#now positions are constraint of the Alanine and the NVT equilibration can be run
	filename = 'nvtInfo'+str(i)+'.dat'
	simulations[i].reporters.append(StateDataReporter(filename, 100, step=True, potentialEnergy=True, temperature=True, separator=' '))
	simulations[i].step(nequib_steps)
print "Equilibration done..."
#Reassigning the masses again
for i in range(Alanine_atoms):
	system.setParticleMass(i,Alanine_masses[i].value_in_unit(dalton))
print "Production started ..."
for iteration in range(max_iterations):
	for i in range(len(simulations)):
		print "We are at temperature " + str(integrators[i].getTemperature())
		simulations[i].step(nsteps_per_iteration)
  		state = simulations[i].context.getState(getEnergy=True, enforcePeriodicBox = True)
    	print "Iteration %5d / %5d | kinetic %8.3f kJ/mol | potential %8.3f kJ/mol" % (iteration, max_iterations, state.getKineticEnergy() / kilojoules_per_mole, state.getPotentialEnergy() / kilojoules_per_mole)
	#mix replicas
	print "We can now exchange the replicas"
	mixReplicas(iteration)


ftime=time.time()

print "total time taken is: "+ str(ftime-stime)










