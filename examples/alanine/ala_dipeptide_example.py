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

import sys, os
sys.path.append(os.path.abspath('../'))
sys.path.append(os.path.abspath('../../'))

 
# in principle, all of these imports should be simplified once this is a
# package
from opentis.orderparameter import OP_Function, OP_Volume
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


from opentis.openmm_simulation import OpenMMSimulation

if __name__=="__main__":
    options = {
                'temperature' : 300.0 * kelvin,
                'collision_rate' : 1.0 / picoseconds,
                'timestep' : 2.0 * femtoseconds,
                'nsteps_per_frame' : 10,
                'n_frames_max' : 5000,
                'start_time' : time.time(),
                'fn_initial_pdb' : "../data/Alanine_solvated.pdb",
                'platform' : 'fastest',
                'solute_indices' : range(22), # TODO: This could be determined automatically !?!?
                'forcefield_solute' : 'amber96.xml',
                'forcefield_solvent' : 'tip3p.xml'
               }
    simulator = OpenMMSimulation(
                    filename="trajectory.nc",
                    topology_file="../data/Alanine_solvated.pdb",
                    opts=options,
                    mode='create'
                    )

    simulator.equilibrate(5)
    snap = Snapshot(simulator)
    simulator.storage.snapshot.save(snap, 0)
    simulator.initialized = True
    PathMover.simulator = simulator

    # this generates an order parameter (callable) object named psi (so if
    # we call `psi(trajectory)` we get a list of the values of psi for each
    # frame in the trajectory). This particular order parameter uses
    # mdtraj's compute_dihedrals function, with the atoms in psi_atoms

    psi_atoms = [6,8,14,16]
    psi = OP_Function("psi", md.compute_dihedrals, trajdatafmt="mdtraj",
                      indices=[psi_atoms])

    # same story for phi, although we won't use that

    phi_atoms = [4,6,8,14]
    phi = OP_Function("phi", md.compute_dihedrals, trajdatafmt="mdtraj",
                      indices=[phi_atoms])

    # save the orderparameters in the storage
    # since they have no data cache this will only contain their name
    psi.save(storage=simulator.storage.cv)
    phi.save(storage=simulator.storage.cv)

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
    for no, interface in enumerate(interface_set):
        # Give each interface a name
        interface.name = 'Interface '+str(no)
        # And save all of these
        simulator.storage.ensemble.save(interface)

    mover_set = mf.OneWayShootingSet(UniformSelector(), interface_set)

    print """
PART ONE: Generate an initial trajectory which satisfies the path ensemble
for the innermost interface.

We do this by using a special sequential ensemble for the sequence.
This path ensemble is particularly complex because we want to be sure that
the path we generate is in the ensemble we desire: this means that we can't
use LeaveXEnsemble as we typically do with TIS paths.
    """
    snapshot = simulator.storage.snapshot.load(0)
    
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
    bootstrap = Bootstrapping(storage=simulator.storage,
                              simulator=simulator,
                              ensembles=interface_set,
                              movers=mover_set)

    bootstrap.set_replicas([segments[0]])

    bootstrap.run(50)

    print """
    Saving all cached computations of orderparameters.
    """

    psi.save(storage=simulator.storage.cv)
    phi.save(storage=simulator.storage.cv)

    # Save all interface volumes as orderparameters
    op_vol_set = [OP_Volume('OP' + str(idx), vol) for idx, vol in enumerate(volume_set)]

    for op in op_vol_set:
        op(simulator.storage.snapshot.all())
        simulator.storage.cv.save(op)

    # Create an orderparameter from a volume
    op_inA = OP_Volume('StateA', stateA)
    op_inB = OP_Volume('StateB', stateB)
    op_notinAorB = OP_Volume('StateX', ~ (stateA | stateB))

    # compute the orderparameter for all snapshots
    op_inA(simulator.storage.snapshot.all())
    op_inB(simulator.storage.snapshot.all())
    op_notinAorB(simulator.storage.snapshot.all())

    simulator.storage.cv.save(op_inA)
    simulator.storage.cv.save(op_inB)
    simulator.storage.cv.save(op_notinAorB)

    # Alternatively one could write

    # simulator.storage.cv.save(psi)
    # simulator.storage.cv.save(phi)
