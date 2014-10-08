"""
@author: David W.H. Swenson
"""
# hack until this is a proper package
import sys
import os
sys.path.append(os.path.abspath('../'))

import numpy as np
 
from Simulator import Simulator
from orderparameter import OP_Function
from snapshot import Snapshot
from volume import LambdaVolume
from ensemble import EnsembleFactory as ef
from ensemble import LengthEnsemble

from simtk.unit import femtoseconds, picoseconds, nanometers, kelvin, dalton
from simtk.unit import Quantity

import time

from integrators import VVVRIntegrator
from simtk.openmm.app.pdbfile import PDBFile
import simtk.openmm as openmm
from simtk.openmm.app import ForceField, PME, HBonds
from simtk.openmm.app import Simulation


# TODO: figure out how much of this should be moved to Simulator 

class AlanineDipeptideTrajectorySimulator(Simulator):
    def __init__(self, filename, topology, opts, mode='auto'):
        opts['fn_storage'] = filename
        self.pdb = PDBFile(topology)
        opts['topology'] = self.pdb.topology

        if mode == 'create':
            platform = openmm.Platform.getPlatformByName(opts['platform'])
            forcefield = ForceField( opts['forcefield_solute'],
                                     opts['forcefield_solvect'] )   
            system = forcefield.createSystem( opts['topology'], 
                                              nonbondedMethod=PME,
                                              nonbondedCutoff=1*nanometers,
                                              constraints=HBonds )
            opts['system_serial'] = openmm.XmlSerializer.serialize(system)
            integrator = VVVRIntegrator( opts['temperature'],
                                         opts['collision_rate'],
                                         opts['timestep'] )
            simulation = Simulation( opts['topology'], system, 
                                     integrator, platform )
            opts['integrator_serial'] = openmm.XmlSerializer.serialize(system)
            opts['integrator_class'] = type(integrator).__name__

            self.max_length_stopper = LengthEnsemble(slice(0,opts['n_frames_max']-1))
            self.simulation = simulation
        if mode == 'restore':
            pass
        self.opts = opts
        pass

    def __getattr__(self, name):
        if name in self.opts.keys():
            return self.opts[name]
        else:
            raise AttributeError

    def equilibrate(self, nsteps):
        self.simulation.context.setPositions(self.pdb.positions)
        system = self.simulation.system
        n_solute = len(self.solute_indices)

        solute_masses = Quantity(np.zeros(n_solute, np.double), dalton)
        for i in self.solute_indices:
            solute_masses[i] = system.getParticleMass(i)
            system.setParticleMass(i,0.0)

        self.simulation.step(nsteps)

        for i in self.solute_indices:
            system.setParticleMass(i, solute_masses[i].value_in_unit(dalton))


    def run_trajectory(self, nsteps):
        pass


if __name__=="__main__":
    start_time = time.time()
    options = {
                'temperature' : 300.0 * kelvin,
                'collision_rate' : 1.0 / picoseconds,
                'timestep' : 2.0 * femtoseconds,
                'nframes_per_iteration' : 10,
                'n_frames_max' : 5000,
                'start_time' : time.time(),
                'fn_initial_pdb' : "../data/Alanine_solvated.pdb",
                'platform' : 'CPU',
                'solute_indices' : range(22),
                'forcefield_solute' : 'amber99sbildn.xml',
                'forcefield_solvect' : 'tip3p.xml'
               }
    simulator = AlanineDipeptideTrajectorySimulator(
                    filename="trajectory.nc",
                    topology="../data/Alanine_solvated.pdb",
                    opts=options,
                    mode='create'
                    )
    
    
    simulator.equilibrate(5)

    snapshot = Snapshot.load(0,0,)
    traj = simulator.generate(snapshot, [LengthEnsemble(slice(0,50))])
    

