"""
This is an example file of what a user might actually need to do to run a
TIS simulation on alanine dipeptide.

@author: David W.H. Swenson
"""
# NOTE: while in development, the goal here is to get something working,
# then to abstract out anything that isn't really necessary. The hope is
# that more of what the user will need to do will be in the __main__ of this
# script

# hack until this is a proper package
import sys
import os
sys.path.append(os.path.abspath('../'))

import numpy as np
import mdtraj as md
 
from Simulator import Simulator
from orderparameter import OP_Function
from snapshot import Snapshot, Configuration
from volume import LambdaVolume
from ensemble import EnsembleFactory as ef
from ensemble import LengthEnsemble
from storage import TrajectoryStorage
from trajectory import Trajectory

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
        super(AlanineDipeptideTrajectorySimulator, self).__init__()
        Snapshot.simulator = self
        Configuration.simulator = self
        Trajectory.simulator = self
        self.opts = {}
        self.add_stored_parameters(opts)
        self.fn_storage = filename 

        self.pdb = PDBFile(topology)
        self.topology = self.pdb.topology

        if mode == 'create':
            self.storage = TrajectoryStorage(
                                    topology=self.topology,
                                    filename=self.fn_storage,
                                    mode='create' )
            platform = openmm.Platform.getPlatformByName(self.platform)
            forcefield = ForceField( self.forcefield_solute,
                                     self.forcefield_solvent )
            system = forcefield.createSystem( self.topology, 
                                              nonbondedMethod=PME,
                                              nonbondedCutoff=1*nanometers,
                                              constraints=HBonds )
            self.system_serial = openmm.XmlSerializer.serialize(system)

            integrator = VVVRIntegrator( self.temperature,
                                         self.collision_rate,
                                         self.timestep )
            self.integrator_serial = openmm.XmlSerializer.serialize(system)
            self.integrator_class = type(integrator).__name__

            simulation = Simulation( self.topology, system, 
                                     integrator, platform )

            self.max_length_stopper = LengthEnsemble(slice(0,self.n_frames_max-1))
            self.simulation = simulation
            self.storage.simulator = self
            self.storage.init_classes()
            self.storage._store_options(self)
        if mode == 'restore':
            pass
        return


    def add_stored_parameters(self, param_dict):
        '''Adds parameters in param_dict to the attribute dictionary of the
        simulator object, and saves the relevant keys as options to store.

        Parameters
        ----------
        param_dict : dict
            dictionary of attributes to be added (and stored); attribute
            names are keys, with appropriate values
        '''
        # TODO: I think this should go into the Simulator object
        for key in param_dict.keys():
            self.opts[key] = param_dict[key]
            setattr(self, key, param_dict[key])
        self.options_to_store = self.opts.keys()
        return
        

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
                'forcefield_solute' : 'amber96.xml',
                'forcefield_solvent' : 'tip3p.xml'
               }
    simulator = AlanineDipeptideTrajectorySimulator(
                    filename="trajectory.nc",
                    topology="../data/Alanine_solvated.pdb",
                    opts=options,
                    mode='create'
                    )
    
    
    simulator.equilibrate(5)
    Snapshot(simulator.simulation.context).save(0,0)
    simulator.initialized = True

    snapshot = Snapshot.load(0,0)
    # this generates an order parameter (callable) object named psi (so if
    # we call `psi(trajectory)` we get a list of the values of psi for each
    # frame in the trajectory). This particular order parameter uses
    # mdtraj's compute_dihedrals function, with the atoms in psi_atoms
    psi_atoms = [6,8,14,16]
    psi = OP_Function("psi", md.compute_dihedrals, trajdatafmt="mdtraj",
                            indices=[psi_atoms])
    # same story for phi, although we won't use that
    phi_atoms = [4,6,8,14]
    phi = OP_Function("psi", md.compute_dihedrals, trajdatafmt="mdtraj",
                            indices=[phi_atoms])

    # now we define our states and our interfaces
    degrees = 180/3.14159 # psi reports in radians; I think in degrees
    stateA = LambdaVolume(psi, -120.0/degrees, -30.0/degrees)
    stateB = LambdaVolume(psi, 100/degrees, 180/degrees) # TODO: periodic?
    interface0 = LambdaVolume(psi, -120.0/degrees, -30.0/degrees)
    
    # TODO: this should be made into some wrapper to clean it up for the
    # user
    paths_AA_interface_0 = ef.TISEnsemble(stateA, stateA, interface0, lazy=True)
    paths_AB_interface_0 = ef.TISEnsemble(stateA, stateB, interface0, lazy=True)

    traj = simulator.generate(snapshot, [LengthEnsemble(slice(0,100))])

    print "Does our full trajectory satisfy the ensemble?",
    print (paths_AA_interface_0(traj) or paths_AB_interface_0(traj))
    for frame in traj:
        print phi(frame)[0], psi(frame)[0], stateA(frame), interface0(frame), stateB(frame)


