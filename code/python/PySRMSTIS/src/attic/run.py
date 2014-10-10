#!/usr/local/bin/env python
#from __main__ import ERROR

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Replica-exchange transition path sampling in the s-field on a Kob-Andersen system.

DESCRIPTION


REFERENCES

[1] Hedges LO, Jack RL, Garrahan JP, and Chandler D. Dynamic order-disorder in atomic models
of structural glass-formers. Science 323:1309, 2009.

[2] Minh DDL and Chodera JD. Optimal estimators and asymptotic variances for nonequilibrium
path-ensemble averages. JCP 131:134110, 2009.

COPYRIGHT

@author John D. Chodera <jchodera@gmail.com>
@author Jan-Hendrik Prinz <jan.prinz@gmx.de>

This source file is released under the GNU General Public License.

This program is free_idx software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import os
import os.path
import sys
import numpy
import math
import copy
import time

from trajectory import Trajectory
from snapshot import Snapshot
from transition_path_sampling import TransitionPathSampling
from replica_exchange_tps import ReplicaExchangeTPS
#import scipy.optimize # necessary for systems with finnicky SciPy installations

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

from simtk.unit import nanosecond, picosecond, nanometers, nanometer, picoseconds, femtoseconds, femtosecond

#import Scientific.IO.NetCDF as netcdf # for netcdf interface in Scientific
import netCDF4 as netcdf # for netcdf interface provided by netCDF4 in enthought python

#=============================================================================================
# SOURCE CONTROL
#=============================================================================================

__version__ = "$Id: KobAndersen.py 524 2010-01-05 07:47:29Z jchodera $"

#=============================================================================================
# GLOBAL CONSTANTS
#=============================================================================================

kB = BOLTZMANN_CONSTANT_kB * AVOGADRO_CONSTANT_NA

#=============================================================================================
# Setup system
#=============================================================================================

def AlanineSetup(N=150, NA=None, A_fraction=0.8, principal_component_density=0.96, mm=None, mass=None, epsilon=None, sigma=None, softcore=False, alpha=0.5, lambda_=1.0):
    """
    Create a test system containing a Kob-Andersen two-component mixture.

    A soft-core Lennard-Jones potential is used if 'softcore' is set to True.

    OPTIONAL ARGUMENTS

    N (int) - total number of atoms (default: 150)
    A_fraction (float) - fraction of A component
    principal_component_density (float) - NA sigma^3 / V (default: 0.96)
    softcore (bool) - soft-core Lennard Jones (Eq. 4 of Shirts and Pande, JCP 122:134508, 2005) is used if True (default: False)
    lambda_ (float) - alchemical parameter, where 1.0 is fully interacting, 0.0 is non-interacting (default: 1.0)
    alpha (float) - soft-core parameter (default: 0.5)

    RETURNS

    system (System)
    coordinates (numpy array)
    epsilon (simtk.unit) - fundamental energy scale (change to argument?)
    sigma (simtk.unit) - fundamental length scale (change to argument?

    EXAMPLES

    Create a Kob-Andersen two-component mixture.

    >>> epsilon     = 119.8 * kelvin * BOLTZMANN_CONSTANT_kB * AVOGADRO_CONSTANT_NA # arbitrary reference energy    
    >>> [system, coordinates] = KobAndersen(epsilon=epsilon)

    Create softcore Kob-Andersen two-component mixture with alchemical perturbation.

    >>> [system, coordinates] = KobAndersen(epsilon=epsilon, softcore=True, lambda_=0.0)

    Test the energy

    >>> # Create a Context.
    >>> kB = BOLTZMANN_CONSTANT_kB
    >>> NA = AVOGADRO_CONSTANT_NA
    >>> temperature = 0.6 * epsilon / kB / NA
    >>> collision_rate = 90.0 / picosecond
    >>> timestep = 1.0 * femtosecond    
    >>> integrator = openLangevinIntegrator(temperature, collision_rate, timestep)
    >>> context = openContext(system, integrator)
    >>> # Set positions
    >>> context.setPositions(coordinates)
    >>> # Evaluate the potential energy.
    >>> state = context.getState(getEnergy=True)
    >>> reduced_potential = (state.getPotentialEnergy() / epsilon)
    >>> print reduced_potential

    Integrate dynamics

    >>> nsteps = 1000 # number of steps to integrate
    >>> integrator.step(nsteps)
    >>> # Retrieve configuration to make sure no coordinates are nan
    >>> state = context.getState(getPositions=True)
    >>> coordinates = state.getPositions(asNumpy=True)
    >>> if numpy.any(numpy.isnan(coordinates / nanometers)): raise Exception('some coordinates are nan after integration: %s' % str(coordinates))

    """

        
    # Create system
    system = System()
    
    principal_component_density = 10.0

    # Compute total system volume.
    volume = NA * sigma**3 / principal_component_density
    
    # Make system cubic in dimension.
    length = volume**(1./3.)
    # TODO: Can we change this to use tuples or 3x3 array?
    a = Quantity(numpy.array([1.0, 0.0, 0.0], numpy.float32), nanometer) * length/nanometer
    b = Quantity(numpy.array([0.0, 1.0, 0.0], numpy.float32), nanometer) * length/nanometer
    c = Quantity(numpy.array([0.0, 0.0, 1.0], numpy.float32), nanometer) * length/nanometer
    print "box edge length = %s" % str(length)
    system.setDefaultPeriodicBoxVectors(a, b, c)
    
    coordinates = []
    
    #
    # TODO: THIS WILL COME FROM BAS
    #
           
    # Return system and coordinates.
    return (system, coordinates)

#=============================================================================================
# CONFIGURATIONAL SPACE SET DEFINITION
#=============================================================================================

# A state must be a decomposition of (potentially overlapping) convex objects in the full configurational space

class ConfigurationalRegion(object):
    """
    ConfigurationalRegion is an object that describes a region in the full configurational space. 
    
    NOTES
    
    It has an associated system object that describes the configurational space. 
    """
    def __init__(self):
        return
    
    def contains(self, state):
        return False
    
    def __call__(self, region):
        return False
    
    
    
#=============================================================================================
# TIS STATE/INTERFACE DEFINITION
#=============================================================================================

class TISState(object):
    """
    CoreState Definition for TIS cores and its associated interfaces
    
    NOTES
    
    Might be reasonable to use these __call__ of a MultiStateDefinition object
    """
    def __init__(self):
        return

#=============================================================================================
# MSM STATE/INTERFACE DEFINITION
#=============================================================================================

class MSMState(object):
    """
    MSM State Definition for TIS cores and its associated interfaces
    
    NOTES
    
    Might be reasonable to use these __call__ of a MultiStateDefinition object
    """
    def __init__(self):
        return

#=============================================================================================
# MAIN AND TESTS
#=============================================================================================

def driver():
    # Constant
    kB = BOLTZMANN_CONSTANT_kB * AVOGADRO_CONSTANT_NA

    # Create a Kob-Andersen two-component mixture.
    N = 150 # number of particles
    A_fraction = 0.8 # fraction of A component
    NA = int(math.floor(A_fraction * N)) # number of A component
    mass        = 39.948 * amu # arbitrary reference mass        
    epsilon     = 119.8 * kelvin * BOLTZMANN_CONSTANT_kB * AVOGADRO_CONSTANT_NA # arbitrary reference energy    
    sigma       = 0.3405 * nanometers # arbitrary reference lengthscale    
    [system, coordinates] = KobAndersen(N=N, NA=NA, mass=mass, epsilon=epsilon, sigma=sigma)    
    
    # Choose the platform and device we will run on.
    #platform = openPlatform.getPlatformByName('Reference') # slow CPU platform
    platform = Platform.getPlatformByName('CPU') # fast GPU platform
    #platform.setPropertyDefaultValue('OpenCLDeviceIndex', '0') # use first OpenCL device

    # Select test or production mode.
    # Test mode runs more quickly, but runs with a different set of parameters.
    #mode = 'test'
    mode = 'production'
    
    # Relevant times and timescales
    reduced_time = sqrt((mass * sigma**2) / (48.0 * epsilon)) # factor on which all reduced times are based
    timestep = 0.035 * reduced_time # velocity Verlet timestep from Ref. [1], Supplementary Section 1.1
    delta_t = (40.0 / 3.0) * reduced_time # number of timesteps per Delta t in Lester's paper, Ref. [1] Supplementary Section 1.1, specified exactly by Lester in private communication

    if mode == 'test':
        quenched_temperature = 0.6 * epsilon / kB # temperature corresponding to left column of Fig. 2 from [1]
        tau = 200.0 * reduced_time # private communication from Lester Hedges (tau = 60 * reduced time for Tred = 0.7, 200 * reduced_time for Tred = 0.6)
        t_obs = 1 * tau # number of intervals of tau, from Fig. 2 of [1] (red dashed line, right column,) # should be 33 * tau for Tred = 0.7, 20*tau for Tred = 0.6
    elif mode == 'production':
        quenched_temperature = 0.7 * epsilon / kB # temperature corresponding to left column of Fig. 2 from [1]
        tau = 60.0 * reduced_time # private communication from Lester Hedges (tau = 60 * reduced time for Tred = 0.7, 200 * reduced_time for Tred = 0.6)
        t_obs = 33 * tau # number of intervals of tau, from Fig. 2 of [1] (red dashed line, right column,) # should be 33 * tau for Tred = 0.7, 20*tau for Tred = 0.6        
    else:
        raise Exception("Unknown mode: %s" % mode)

    print "UNCORRECTED TIMES"
    print "timestep = %s" % str(timestep.in_units_of(femtosecond))
    print "delta_t = %s (%f steps)" % (str(delta_t.in_units_of(picosecond)), delta_t / timestep)
    print "tau = %s (%f steps)" % (str(tau.in_units_of(picosecond)), tau / timestep)
    print "t_obs = %s (%f steps)" % (str(t_obs.in_units_of(picosecond)), t_obs / timestep)

    # Determine integral number of steps per velocity randomization and number of trajectories
    nsteps_per_frame = int(round(delta_t / timestep)) # number of steps per trajectory frame and velocity randomization
    nframes = int(round(t_obs / delta_t)) # number of frames per trajectory
    print "number of steps per trajectory frame and velocity randomization = %d" % nsteps_per_frame
    print "number of frames per trajectory = %d" % nframes
    # Correct delta_t and t_obs to be integral
    delta_t = nsteps_per_frame * timestep 
    t_obs = nframes * delta_t
    print "CORRECTED TIMES (after rounding number of steps to integers)"
    print "timestep = %s" % str(timestep.in_units_of(femtosecond))
    print "delta_t = %s (%f steps)" % (str(delta_t.in_units_of(picosecond)), delta_t / timestep)
    print "tau = %s (%f steps)" % (str(tau.in_units_of(picosecond)), tau / timestep)
    print "t_obs = %s (%f steps)" % (str(t_obs.in_units_of(picosecond)), t_obs / timestep)

    # Replica-exchange filename
    ncfilename = 'repex.nc'

    # Set default simulation protocol.
    minimize = True
    equilibrate = True
    quench = True
    seed = True
    
    # If we're resuming from a previously-existing NetCDF file, change the protocol to skip equilibration and quenching.
    if os.path.exists(ncfilename):
        # No need to do these things if we are resuming
        minimize = True
        equilibrate = False
        quench = False
        seed = True
                
    if minimize:
        # Minimize the system prior to dynamics.
        print "Minimizing with L-BFGS..."
        # Initialize a minimizer with default options.
        minimizer = LocalEnergyMinimizer(system, verbose=True, platform=platform)
        # Minimize the initial coordinates.
        
        coordinates = minimizer.minimize(coordinates)
        # Clean up to release the Context.
        del minimizer


    # Set temperature for equilibration simulations.
    elevated_temperature = 2.0 * epsilon / kB

    if equilibrate:
        # Equilibrate at high temperature.
        collision_rate = 1.0 / delta_t
        nsteps = int(math.floor(t_obs / timestep))
        print "Equilibrating at %s for %d steps..." % (str(elevated_temperature), nsteps)
        integrator = LangevinIntegrator(elevated_temperature, collision_rate, timestep)
        context = Context(system, integrator, platform)
        context.setPositions(coordinates)
        integrator.step(nsteps)
        state = context.getState(getPositions=True)
        coordinates = state.getPositions(asNumpy=True)
        del state, context, integrator

    if quench:
        # Quench to final temperature
        collision_rate = 1.0 / delta_t
        nsteps = int(math.floor(t_obs / timestep))
        print "Quenching to %s for %d steps..." % (str(quenched_temperature), nsteps)
        integrator = LangevinIntegrator(quenched_temperature, collision_rate, timestep)
        context = Context(system, integrator, platform)
        context.setPositions(coordinates)
        integrator.step(nsteps)
        state = context.getState(getPositions=True)
        coordinates = state.getPositions(asNumpy=True)
        del state, context, integrator

    # Create a number of transition path sampling ensembles at different values of the field parameter s.
    s_reduced_unit = 1.0 / (sigma**2 * delta_t)
    #svalues = [0.00, 0.01, 0.02, 0.03, 0.04, 0.06] # minimal set of s values
    svalues = [0.00, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, -0.005]
    print "Initializing transition path sampling ensembles..."
    ensembles = [ TransitionPathSampling(system, N, NA, mass, epsilon, sigma, timestep, nsteps_per_frame, nframes, quenched_temperature, s * s_reduced_unit, platform=platform) for s in svalues ]

    trajectory = None
    if seed:
        # Generate an initial trajectory at zero field
        print "Generating seed trajectory for TPS..."
        trajectory = ensembles[0].generateTrajectory(coordinates, nframes)

    # Initialize replica-exchange TPS simulation.
    print "Initializing replica-exchange TPS..."
    simulation = ReplicaExchangeTPS(system, ensembles, trajectory, ncfilename)

    # Set simulation parameters.
    simulation.number_of_iterations = 10000

    # Run simulation
    print "Running replica-exchange TPS..."
    simulation.run()
        
    print "Work Done :) Go home and have a beer."

#=============================================================================================
# This driver sets up a simulation for multiple temperatures and s-values.
#=============================================================================================
def multitemp_driver():
    # Constant
    kB = BOLTZMANN_CONSTANT_kB * AVOGADRO_CONSTANT_NA

    # Create a Kob-Andersen two-component mixture.
    N = 150 # number of particles
    A_fraction = 0.8 # fraction of A component
    NA = int(math.floor(A_fraction * N)) # number of A component
    mass        = 39.948 * amu # arbitrary reference mass        
    epsilon     = 119.8 * kelvin * BOLTZMANN_CONSTANT_kB * AVOGADRO_CONSTANT_NA # arbitrary reference energy    
    sigma       = 0.3405 * nanometers # arbitrary reference lengthscale    
    [system, coordinates] = KobAndersen(N=N, NA=NA, mass=mass, epsilon=epsilon, sigma=sigma)    

    # Choose the platform and device we will run on.
    #platform = openPlatform.getPlatformByName('Reference') # slow CPU platform
    platform = Platform.getPlatformByName('CPU') # fast GPU platform
    #platform.setPropertyDefaultValue('OpenCLDeviceIndex', '0') # use first OpenCL device

    # Select test or production mode.
    # Test mode runs more quickly, but runs with a different set of parameters.
    #mode = 'test'
    mode = 'production'
    
    # Relevant times and timescales
    reduced_time = sqrt((mass * sigma**2) / (48.0 * epsilon)) # factor on which all reduced times are based
    timestep = 0.035 * reduced_time # velocity Verlet timestep from Ref. [1], Supplementary Section 1.1
    delta_t = (40.0 / 3.0) * reduced_time # number of timesteps per Delta t in Lester's paper, Ref. [1] Supplementary Section 1.1, specified exactly by Lester in private communication

    if mode == 'test':
        quenched_temperature = 0.6 * epsilon / kB # temperature corresponding to left column of Fig. 2 from [1]
        tau = 200.0 * reduced_time # private communication from Lester Hedges (tau = 60 * reduced time for Tred = 0.7, 200 * reduced_time for Tred = 0.6)
        t_obs = 1 * tau # number of intervals of tau, from Fig. 2 of [1] (red dashed line, right column,) # should be 33 * tau for Tred = 0.7, 20*tau for Tred = 0.6
    elif mode == 'production':
        quenched_temperature = 0.7 * epsilon / kB # temperature corresponding to left column of Fig. 2 from [1]
        tau = 60.0 * reduced_time # private communication from Lester Hedges (tau = 60 * reduced time for Tred = 0.7, 200 * reduced_time for Tred = 0.6)
        t_obs = 33 * tau # number of intervals of tau, from Fig. 2 of [1] (red dashed line, right column,) # should be 33 * tau for Tred = 0.7, 20*tau for Tred = 0.6        
    else:
        raise Exception("Unknown mode: %s" % mode)

    print "UNCORRECTED TIMES"
    print "timestep = %s" % str(timestep.in_units_of(femtosecond))
    print "delta_t = %s (%f steps)" % (str(delta_t.in_units_of(picosecond)), delta_t / timestep)
    print "tau = %s (%f steps)" % (str(tau.in_units_of(picosecond)), tau / timestep)
    print "t_obs = %s (%f steps)" % (str(t_obs.in_units_of(picosecond)), t_obs / timestep)

    # Determine integral number of steps per velocity randomization and number of trajectories
    nsteps_per_frame = int(round(delta_t / timestep)) # number of steps per trajectory frame and velocity randomization
    nframes = int(round(t_obs / delta_t)) # number of frames per trajectory
    print "number of steps per trajectory frame and velocity randomization = %d" % nsteps_per_frame
    print "number of frames per trajectory = %d" % nframes
    # Correct delta_t and t_obs to be integral
    delta_t = nsteps_per_frame * timestep 
    t_obs = nframes * delta_t
    print "CORRECTED TIMES (after rounding number of steps to integers)"
    print "timestep = %s" % str(timestep.in_units_of(femtosecond))
    print "delta_t = %s (%f steps)" % (str(delta_t.in_units_of(picosecond)), delta_t / timestep)
    print "tau = %s (%f steps)" % (str(tau.in_units_of(picosecond)), tau / timestep)
    print "t_obs = %s (%f steps)" % (str(t_obs.in_units_of(picosecond)), t_obs / timestep)

    # Replica-exchange filename
    ncfilename = 'repex-Ts.nc' # DEBUG

    # Set default simulation protocol.
    minimize = False
    equilibrate = True
    quench = True
    seed = True
    
    # If we're resuming from a previously-existing NetCDF file, change the protocol to skip equilibration and quenching.
    if os.path.exists(ncfilename):
        # No need to do these things if we are resuming
        minimize = True
        equilibrate = False
        quench = False
        seed = True
                
    if minimize:
        # Minimize the system prior to dynamics.
        print "Minimizing with L-BFGS..."
        # Initialize a minimizer with default options.
        minimizer = optimize.LBFGSMinimizer(system, verbose=True, platform=platform)
        # Minimize the initial coordinates.
        coordinates = minimizer.minimize(coordinates)
        # Clean up to release the Context.
        del minimizer

    # Set temperature for equilibration simulations.
    elevated_temperature = 2.0 * epsilon / kB

    if equilibrate:
        # Equilibrate at high temperature.
        collision_rate = 1.0 / delta_t
        nsteps = int(math.floor(t_obs / timestep))
        print "Equilibrating at %s for %d steps..." % (str(elevated_temperature), nsteps)
        integrator = LangevinIntegrator(elevated_temperature, collision_rate, timestep)
        context = Context(system, integrator, platform)
        context.setPositions(coordinates)
        integrator.step(nsteps)
        state = context.getState(getPositions=True)
        coordinates = state.getPositions(asNumpy=True)
        del state, context, integrator

    if quench:
        # Quench to final temperature
        collision_rate = 1.0 / delta_t
        nsteps = int(math.floor(t_obs / timestep))
        print "Quenching to %s for %d steps..." % (str(quenched_temperature), nsteps)
        integrator = LangevinIntegrator(quenched_temperature, collision_rate, timestep)
        context = Context(system, integrator, platform)
        context.setPositions(coordinates)
        integrator.step(nsteps)
        state = context.getState(getPositions=True)
        coordinates = state.getPositions(asNumpy=True)
        del state, context, integrator

    # Create a number of transition path sampling ensembles at different values of the field parameter s and different temperatures
    s_reduced_unit = 1.0 / (sigma**2 * delta_t) # reduced units for s field
    T_reduced_unit = epsilon / kB # reduced units for temperature

    svalues = [0.00, 0.01, 0.02, 0.03, 0.04, 0.06] # minimal set of s values
    Tvalues = [0.70, 0.80, 0.90, 1.00, 1.10, 1.20] # some made-up temperatures
               
    print "Initializing transition path sampling ensembles..."
    ensembles = [ TransitionPathSampling(system, N, NA, mass, epsilon, sigma, timestep, nsteps_per_frame, nframes, T * T_reduced_unit, s * s_reduced_unit, platform=platform) for (s,T) in zip(svalues,Tvalues) ]

    trajectory = None
    if seed:
        # Generate an initial trajectory at zero field
        print "Generating seed trajectory for TPS..."
        trajectory = ensembles[0].generateTrajectory(coordinates, nframes)

    # Initialize replica-exchange TPS simulation.
    print "Initializing replica-exchange TPS..."
    simulation = ReplicaExchangeTPS(system, ensembles, trajectory, ncfilename)

    # Set simulation parameters.
    simulation.number_of_iterations = 10

    # Run simulation
    print "Running replica-exchange TPS..."
    print "(I'll buy you a beer if you actually get this to work.)"
    simulation.run()

    print "And there was much rejoicing!"

if __name__ == "__main__":
    #import doctest
    #doctest.testmod()

    # Run replica-exchange TPS driver
    driver()

    # Run new-and-improved multiple temperature and s value driver.
    #multitemp_driver()

    
   
  
