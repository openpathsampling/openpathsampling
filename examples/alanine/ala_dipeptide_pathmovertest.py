"""
This is an example file of what a user might actually need to do to run a
TIS simulation on alanine dipeptide.

@author: David W.H. Swenson
"""
# NOTE: while in development, the goal here is to get something working,
# then to abstract out anything that isn't really necessary. The hope is
# that more of what the user will need to do will be in the __main__ of this
# script

import logging.config
import numpy as np
import mdtraj as md
import time

import sys, os
sys.path.append(os.path.abspath('../'))
sys.path.append(os.path.abspath('../../'))

# in principle, all of these imports should be simplified once this is a
# package
from openpathsampling.orderparameter import OP_Function, OP_Volume, OP_MD_Function
from openpathsampling.openmm_engine import OpenMMEngine
from openpathsampling.snapshot import Snapshot
from openpathsampling.volume import LambdaVolumePeriodic, VolumeFactory as vf
from openpathsampling.pathmover import PathMoverFactory as mf
from openpathsampling.ensemble import EnsembleFactory as ef
from openpathsampling.ensemble import (LengthEnsemble, SequentialEnsemble, OutXEnsemble,
                              InXEnsemble)
from openpathsampling.calculation import Bootstrapping
from openpathsampling.pathmover import PathMover, MoveDetails, SequentialMover, \
    ConditionalSequentialMover, PartialAcceptanceSequentialMover, \
    ForwardShootMover, CollapseMove
from openpathsampling.shooting import UniformSelector
from openpathsampling.sample import Sample, SampleSet

import simtk.unit as u

if __name__=="__main__":
    logging.config.fileConfig('../../openpathsampling/logging.conf',
                              disable_existing_loggers=False)
    options = {'temperature' : 300.0 * u.kelvin,
               'collision_rate' : 1.0 / u.picoseconds,
               'timestep' : 2.0 * u.femtoseconds,
               'nsteps_per_frame' : 10,
               'n_frames_max' : 5000,
               'start_time' : time.time(),
               'fn_initial_pdb' : "../data/Alanine_solvated.pdb",
               'platform' : 'fastest',
               'solute_indices' : range(22), # TODO: This could be determined automatically !?!?
               'forcefield_solute' : 'amber96.xml',
               'forcefield_solvent' : 'tip3p.xml'
              }

    engine = OpenMMEngine.auto(
        filename="trajectory.nc",
        template='../data/Alanine_solvated.pdb',
        options=options,
        mode='create'
    )

    # set up the initial conditions
    init_pdb = md.load(options['fn_initial_pdb'], frame=0)
    init_pos = init_pdb.xyz[0]
    init_box = init_pdb.unitcell_vectors[0]
    init_vel = np.zeros(init_pos.shape) # NOTE: need ways to assign Boltzmann
    engine.current_snapshot = Snapshot(coordinates=init_pos,
                                       box_vectors=init_box,
                                       velocities=init_vel)

    engine.equilibrate(5)
    snap = engine.current_snapshot
    engine.storage.snapshot.save(snap, 0)
    engine.initialized = True
    PathMover.engine = engine

    # this generates an order parameter (callable) object named psi (so if
    # we call `psi(trajectory)` we get a list of the values of psi for each
    # frame in the trajectory). This particular order parameter uses
    # mdtraj's compute_dihedrals function, with the atoms in psi_atoms

    psi_atoms = [6,8,14,16]
    psi = OP_MD_Function("psi", md.compute_dihedrals, indices=[psi_atoms])

    # same story for phi, although we won't use that

    phi_atoms = [4,6,8,14]
    phi = OP_MD_Function("phi", md.compute_dihedrals, indices=[phi_atoms])

    # save the orderparameters in the storage
    # since they have no data cache this will only contain their name
    psi.save(storage=engine.storage.collectivevariable)
    phi.save(storage=engine.storage.collectivevariable)

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
        engine.storage.ensemble.save(interface)

    mover_set = mf.OneWayShootingSet(UniformSelector(), interface_set)

    print """
PART ONE: Generate an initial trajectory which satisfies the path ensemble
for the innermost interface.

We do this by using a special sequential ensemble for the sequence.
This path ensemble is particularly complex because we want to be sure that
the path we generate is in the ensemble we desire: this means that we can't
use LeaveXEnsemble as we typically do with TIS paths.
    """
    snapshot = engine.storage.snapshot.load(0)

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

    total_path = engine.generate(snapshot,
                                    [first_traj_ensemble.can_append])

    segments = interface0_ensemble.split(total_path)

    first_sample = Sample(
        replica=0,
        trajectory=segments[0],
        ensemble=interface0_ensemble,
        details=MoveDetails()
    )

    first_set = SampleSet([first_sample])

    print first_set[0].__dict__
    print 'Set #1', first_set.__dict__
    print

    mover = ForwardShootMover(UniformSelector(), ensembles=interface0_ensemble)

    first_path = mover_set[0].move(first_set)

    print first_path.__dict__
    print first_path

    print 'Collapsed'

    print first_path.collapsed_samples


    second_set = first_set + first_path
#    second_set = first_set + first_path

    print second_set.__dict__

    print second_set.__class__

    seq_mover = SequentialMover([mover_set[0]] * 2)
    second_path = seq_mover.move(second_set)

    third_set = second_set + second_path

    print second_path.accepted

    print interface0_ensemble(third_set[0].trajectory)

    print third_set[0].details.__dict__
    print third_set[0].__dict__
    print third_set.__dict__

    print second_path.samples
    print first_path.samples


    mover3 = CollapseMove(SequentialMover([
        mover,
        PartialAcceptanceSequentialMover([mover] * 2),
        CollapseMove(ConditionalSequentialMover([
            PartialAcceptanceSequentialMover([mover] * 3)] * 2
        )),
        mover
    ] * 3))

    third_path = mover3.move(third_set)

    print third_path.samples
    print str(third_path)
    print str(third_path.opened)
    print len(third_path)

    forth_set = third_set + third_path
    print interface0_ensemble(forth_set[0].trajectory)

    engine.storage.sampleset.save(forth_set)

    exit()

    bootstrap = Bootstrapping(storage=engine.storage,
                              engine=engine,
                              ensembles=interface_set,
                              movers=mover_set,
                              trajectory=segments[0]
                             )



    bootstrap.run(50)

    print """
    Saving all cached computations of orderparameters.
    """

    engine.storage.collectivevariable.sync(psi)
    engine.storage.collectivevariable.sync(phi)

    # Save all interface volumes as orderparameters
    op_vol_set = [OP_Volume('OP' + str(idx), vol) for idx, vol in enumerate(volume_set)]

    for op in op_vol_set:
        op(engine.storage.snapshot.all())
        engine.storage.collectivevariable.save(op)

    # Create an orderparameter from a volume
    op_inA = OP_Volume('StateA', stateA)
    op_inB = OP_Volume('StateB', stateB)
    op_notinAorB = OP_Volume('StateX', ~ (stateA | stateB))

    # compute the orderparameter for all snapshots
    op_inA(engine.storage.snapshot.all())
    op_inB(engine.storage.snapshot.all())
    op_notinAorB(engine.storage.snapshot.all())

    engine.storage.collectivevariable.save(op_inA)
    engine.storage.collectivevariable.save(op_inB)
    engine.storage.collectivevariable.save(op_notinAorB)