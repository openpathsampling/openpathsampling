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
#import sys
#import os
#sys.path.append(os.path.abspath('../'))

import numpy as np
import mdtraj as md
 
# in principle, all of these imports should be simplified once this is a
# package
from opentis.Simulator import Simulator
from opentis.orderparameter import OP_Function
from opentis.openmm_simulation import OpenMMSimulation
from opentis.snapshot import Snapshot, Configuration
from opentis.volume import LambdaVolumePeriodic, VolumeFactory as vf
from opentis.pathmover import PathMoverFactory as mf
from opentis.ensemble import EnsembleFactory as ef
from opentis.ensemble import (LengthEnsemble, SequentialEnsemble, OutXEnsemble,
                      InXEnsemble)
from opentis.storage import Storage
from opentis.trajectory import Trajectory
from opentis.calculation import Bootstrapping
from opentis.pathmover import (PathMover, MixedMover, ForwardShootMover, 
                       BackwardShootMover)
from opentis.shooting import UniformSelector

from simtk.unit import femtoseconds, picoseconds, nanometers, kelvin, dalton
from simtk.unit import Quantity

import time

from opentis.integrators import VVVRIntegrator

from simtk.openmm.app.pdbfile import PDBFile
import simtk.openmm as openmm
from simtk.openmm.app import ForceField, PME, HBonds
from simtk.openmm.app import Simulation


# TODO: figure out how much of this should be moved to Simulator 

class AlanineDipeptideTrajectorySimulator(Simulator):
    def __init__(self, filename, topology_file, opts, mode='auto'):
        super(AlanineDipeptideTrajectorySimulator, self).__init__()

        # tell everybody who their simulator is
        Snapshot.simulator = self
        Configuration.simulator = self
        Trajectory.simulator = self

        # set up the opts
        self.opts = {}
        self.add_stored_parameters(opts)

        # storage
        self.fn_storage = filename 

        # topology
        self.pdb = PDBFile(topology_file)
        self.topology = self.pdb.topology

        if mode == 'create':
            # set up the OpenMM simulation
            #self.platform = 'CUDA'
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

            simulation = OpenMMSimulation(self.topology, system, 
                                          integrator, platform )

            # claim the OpenMM simulation as our own
            self.simulation = simulation

            # set up the max_length_stopper (if n_frames_max is given)
            self.max_length_stopper = LengthEnsemble(slice(0,self.n_frames_max-1))

            # storage
            self.storage = Storage(
                topology_file=self.topology,
                filename=self.fn_storage,
                mode='w'
            )
            self.storage.simulator = self
            self.storage._store_options(self)
            Trajectory.storage = self.storage
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
        # TODO: rename... this is position restrained equil, right?
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
    options = {
                'temperature' : 300.0 * kelvin,
                'collision_rate' : 1.0 / picoseconds,
                'timestep' : 2.0 * femtoseconds,
                'nframes_per_iteration' : 50,
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
                    topology_file="../data/Alanine_solvated.pdb",
                    opts=options,
                    mode='create'
                    )
    
    
    simulator.equilibrate(5)
    snap = Snapshot(simulator.simulation)
    simulator.storage.snapshot.save(snap, 0, 0)
    simulator.initialized = True
    PathMover.simulator = simulator

    # this generates an order parameter (callable) object named psi (so if
    # we call `psi(trajectory)` we get a list of the values of psi for each
    # frame in the trajectory). This particular order parameter uses
    # mdtraj's compute_dihedrals function, with the atoms in psi_atoms
    psi_atoms = [6,8,14,16]
    psi = OP_Function("psi", md.compute_dihedrals, trajdatafmt="mdtraj",
                      indices=[psi_atoms],
                      storages=simulator.storage.configuration)
    # same story for phi, although we won't use that
    phi_atoms = [4,6,8,14]
    phi = OP_Function("phi", md.compute_dihedrals, trajdatafmt="mdtraj",
                      indices=[phi_atoms],
                      storages=simulator.storage.configuration)

    psi.save()
    phi.save()

    # now we define our states and our interfaces
    degrees = 180/3.14159 # psi reports in radians; I think in degrees
    stateA = LambdaVolumePeriodic(psi, -120.0/degrees, -30.0/degrees)
    stateB = LambdaVolumePeriodic(psi, 100/degrees, 180/degrees) 

    # set up minima and maxima for this transition's interface set
    minima = map((1.0 / degrees).__mul__,
                 [-125, -135, -140, -142.5, -145.0, -147.0, 150.0])
    maxima = map((1.0 / degrees).__mul__,
                 [-25.0, -21.0, -18.5, -17.0, -15.0, -10.0, 0.0])

    volume_set = vf.LambdaVolumePeriodicSet(psi, minima, maxima)
    interface0 = volume_set[0]
    interface_set = ef.TISEnsembleSet(stateA, stateA | stateB, volume_set)
    mover_set = mf.OneWayShootingSet(UniformSelector(), interface_set)

    print """
PART ONE: Generate an initial trajectory which satisfies the path ensemble
for the innermost interface.

We do this by using a special sequential ensemble for the sequence.
This path ensemble is particularly complex because we want to be sure that
the path we generate is in the ensemble we desire: this means that we can't
use LeaveXEnsemble as we typically do with TIS paths.
    """
    snapshot = simulator.storage.snapshot.load(0,0)

    
    first_traj_ensemble = SequentialEnsemble([
        OutXEnsemble(stateA) | LengthEnsemble(0),
        InXEnsemble(stateA),
        (OutXEnsemble(stateA) & InXEnsemble(interface0)) | LengthEnsemble(0),
        InXEnsemble(interface0) | LengthEnsemble(0),
        OutXEnsemble(interface0),
        OutXEnsemble(stateA) | LengthEnsemble(0),
        InXEnsemble(stateA) & LengthEnsemble(1)
    ]) 

    interface0_ensemble = interface_set[0]
    print "start path generation (should not take more than a few minutes)"
    total_path = simulator.generate(snapshot, [first_traj_ensemble.forward])
    print "path generation complete"
    print
    print "Total trajectory length: ", len(total_path)
    segments = interface0_ensemble.split(total_path)
    print "Traj in first_traj_ensemble? (should be)", 
    print first_traj_ensemble(total_path)
    print "Traj in TIS ensemble? (probably not)", 
    print interface0_ensemble(total_path)
    print "Number of segments in TIS ensemble: ", len(segments)
    if len(segments):
        print "Length of each segment:"
        for i in range(len(segments)):
            print "  seg[{0}]: {1}".format(i, len(segments[i]))

    print "Full traj: phi psi ",
    print "stateA iface0 stateB ",
    print "can_append "
    for frame in total_path:
        print phi(frame)[0]*degrees, psi(frame)[0]*degrees, 
        print stateA(frame), interface0(frame), stateB(frame),
        print first_traj_ensemble.forward(
            total_path[slice(0,total_path.index(frame)+1)]
        )

    print
    if len(segments):
        print "Do our segments satisfy the ensemble?",
        for seg in segments:
            print interface0_ensemble(seg),
        print
        print "Segment trajectory: phi psi stateA interface0 stateB"
        for frame in segments[0]:
            print phi(frame)[0]*degrees, psi(frame)[0]*degrees, stateA(frame), interface0(frame), stateB(frame)

    print """
Starting the bootstrapping procedure to obtain initial paths. First we
define our shooting movers (randomly pick fwd or bkwd shooting), then build
the bootstrapping calculation, then we run it. 
    """
    bootstrap = Bootstrapping(simulator.storage, simulator, interface_set,
                              mover_set)
    bootstrap.replicas = [segments[0]]

    bootstrap.run(20)
