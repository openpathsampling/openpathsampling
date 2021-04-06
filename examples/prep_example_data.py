"""
Script to prepare datafiles for analysis.

Occasionally, we need to create a new version of the datafiles that are used
in our continuous integration testing. Typically, this occurs when a new
version of Python changing something about the way bytecode is stored,
meaning that old CVs can't be reloaded in the new Python version.

This script generates the files for the toy MSTIS and toy MISTIS examples.
The resulting files have 1000 MC steps and start with all interfaces
populated by a trajectory that crosses the outermost interface.

Usage:
    python prep_example_data.py

It will determine your Python version and make files named
``toy_{example}_1k_OPS1_py{pyVersion}.nc``, where ``example`` is mstis or
mistis, and `pyVersion` is the major/minor version of Python used (e.g., 39
for Python 3.9).
"""

import itertools
import sys
import numpy as np
import argparse
import openpathsampling as paths
import openpathsampling.engines.toy as toys


def make_mstis_engine():
    pes = (
        toys.OuterWalls([1.0, 1.0], [0.0, 0.0])
        + toys.Gaussian(-0.7, [12.0, 12.0], [0.0, 0.4])
        + toys.Gaussian(-0.7, [12.0, 12.0], [-0.5, -0.5])
        + toys.Gaussian(-0.7, [12.0, 12.0], [0.5, -0.5])
    )

    topology=toys.Topology(
        n_spatial=2,
        masses=[1.0, 1.0],
        pes=pes
    )

    integ = toys.LangevinBAOABIntegrator(dt=0.02, temperature=0.1, gamma=2.5)

    options={
        'integ': integ,
        'n_frames_max': 5000,
        'n_steps_per_frame': 1
    }

    toy_eng = toys.Engine(
        options=options,
        topology=topology
    ).named('toy_engine')

    template = toys.Snapshot(
        coordinates=np.array([[-0.5, -0.5]]),
        velocities=np.array([[0.0,0.0]]),
        engine=toy_eng
    )

    toy_eng.current_snapshot = template
    return toy_eng


def circle(snapshot, center):
    import math
    return math.sqrt((snapshot.xyz[0][0]-center[0])**2
                     + (snapshot.xyz[0][1]-center[1])**2)


def make_mstis_network():
    opA = paths.CoordinateFunctionCV(name="opA", f=circle, center=[-0.5, -0.5])
    opB = paths.CoordinateFunctionCV(name="opB", f=circle, center=[0.5, -0.5])
    opC = paths.CoordinateFunctionCV(name="opC", f=circle, center=[0.0, 0.4])
    stateA = paths.CVDefinedVolume(opA, 0.0, 0.2).named("A")
    stateB = paths.CVDefinedVolume(opB, 0.0, 0.2).named("B")
    stateC = paths.CVDefinedVolume(opC, 0.0, 0.2).named("C")

    interfacesA = paths.VolumeInterfaceSet(opA, 0.0, [0.2, 0.3, 0.4])
    interfacesB = paths.VolumeInterfaceSet(opB, 0.0, [0.2, 0.3, 0.4])
    interfacesC = paths.VolumeInterfaceSet(opC, 0.0, [0.2, 0.3, 0.4])

    ms_outers = paths.MSOuterTISInterface.from_lambdas(
        {ifaces: 0.5
         for ifaces in [interfacesA, interfacesB, interfacesC]}
    )

    mstis = paths.MSTISNetwork(
        [(stateA, interfacesA),
         (stateB, interfacesB),
         (stateC, interfacesC)],
        ms_outers=ms_outers
    )
    return mstis

def mstis_setup():
    engine = make_mstis_engine()
    network = make_mstis_network()
    snapshots = [
        toys.Snapshot(coordinates=np.array([[-0.5, -0.5]]),
                      velocities=np.array([[1.0, 0.0]])),
        toys.Snapshot(coordinates=np.array([[0.5, -0.5]]),
                      velocities=np.array([[-1.0, 0.0]])),
        toys.Snapshot(coordinates=np.array([[0.0, 0.4]]),
                      velocities=np.array([[0.0, -0.5]])),
    ]
    return engine, network, snapshots


def make_mistis_engine():
    pes = (
        toys.OuterWalls([1.0, 1.0], [0.0, 0.0])
        + toys.Gaussian(-1.0, [12.0, 12.0], [-0.5, 0.5])
        + toys.Gaussian(-1.0, [12.0, 12.0], [-0.5, -0.5])
        + toys.Gaussian(-1.0, [12.0, 12.0], [0.5, -0.5])
    )

    topology=toys.Topology(n_spatial=2,
                           masses=[1.0, 1.0],
                           pes=pes)

    integ = toys.LangevinBAOABIntegrator(dt=0.02, temperature=0.1, gamma=2.5)

    options = {
        'integ': integ,
        'n_frames_max': 5000,
        'n_steps_per_frame': 1
    }

    toy_eng = toys.Engine(
        options=options,
        topology=topology
    ).named('engine')
    return toy_eng

def xval(snapshot):
    return snapshot.xyz[0][0]

def xprime(snapshot):
    # this only exists until we set up the ability for the order parameter to decrease
    return -snapshot.xyz[0][0]

def yval(snapshot):
    return snapshot.xyz[0][1]

def make_mistis_network():
    cvX = paths.FunctionCV(name="cvX", f=xval)
    cvY = paths.FunctionCV(name="cvY", f=yval)
    cvXprime = paths.FunctionCV(name="cvXprime", f=xprime)

    x_under_min = paths.CVDefinedVolume(cvX, float("-inf"), -0.35)
    x_over_max = paths.CVDefinedVolume(cvX, 0.35, float("inf"))
    y_under_min = paths.CVDefinedVolume(cvY, float("-inf"), -0.35)
    y_over_max = paths.CVDefinedVolume(cvY, 0.35, float("inf"))

    stateA = (x_under_min & y_under_min).named("A")
    stateB = (x_over_max & y_under_min).named("B")
    stateC = (x_under_min & y_over_max).named("C")

    interfacesAB = paths.VolumeInterfaceSet(
        cvX, float("-inf"), [-0.35, -0.3, -0.27, -0.24, -0.2, -0.1]
    )
    interfacesAC = paths.VolumeInterfaceSet(
        cvY, float("-inf"), [-0.35, -0.3, -0.27, -0.24, -0.2, -0.1, 0.0]
    )
    interfacesBA = paths.VolumeInterfaceSet(
        cvXprime, float("-inf"), [-0.35, -0.3, -0.27, -0.24, -0.2, -0.1]
    )

    ms_outer = paths.MSOuterTISInterface.from_lambdas(
        {iface: 0.0 for iface in [interfacesAB, interfacesBA]}
    )
    network = paths.MISTISNetwork(
        [(stateA, interfacesAB, stateB),
         (stateA, interfacesAC, stateC),
         (stateB, interfacesBA, stateA)],
        ms_outers=ms_outer,
        strict_sampling=True
    ).named("mistis")
    return network

def mistis_setup():
    engine = make_mistis_engine()
    network = make_mistis_network()
    snapshots = [
        toys.Snapshot(coordinates=np.array([[-0.5, -0.5]]),
                      velocities=np.array([[0.5, 0.0]])),
        toys.Snapshot(coordinates=np.array([[0.5, -0.5]]),
                      velocities=np.array([[-0.5, 0.0]]))
    ]
    return engine, network, snapshots

def shoot_until_A_to_A(initial_ensemble, desired_ensemble, sample, engine):
    # we only shoot forward because we know the final frame is the problem
    mover = paths.ForwardShootMover(ensemble=initial_ensemble,
                                    selector=paths.UniformSelector(),
                                    engine=engine)
    while not desired_ensemble(sample):
        change = mover.move_core([sample])
        if desired_ensemble(change.trials[0]):
            sample = change.trials[0]

    return sample

def fill_minus(minus_ensemble, innermost_ensemble, forbidden_states,
               initial_trj, engine):
    initial_state = minus_ensemble.state_vol
    print(f"Filling minus ensemble for state {initial_state.name}")
    forbidden_states_ensemble = paths.AllOutXEnsemble(
        paths.join_volumes(forbidden_states)
    )
    desired_ensemble = innermost_ensemble & forbidden_states_ensemble
    # ensure we're A->A, not A->B
    sample_A_to_A = shoot_until_A_to_A(innermost_ensemble, desired_ensemble,
                                       initial_trj, engine)

    # with an A->A segment, just use this to extend into the minus ensemble
    sample = minus_ensemble.extend_sample_from_trajectories(
        sample_A_to_A,
        engine=engine,
        replica=-1
    )
    return sample.trajectory


def bootstrap(sampling_transition, init_frame, forbidden_states, engine,
              ms_outers):
    my_state = sampling_transition.stateA
    print("Ratcheting initial conditions for state " + my_state.name)
    my_interfaces = sampling_transition.interfaces
    try:
        outer_interface = ms_outers.volume_for_interface_set(my_interfaces)
    except KeyError:
        extra_interfaces = None
    else:
        extra_interfaces = [outer_interface]

    ratchet = paths.FullBootstrapping(
        transition=sampling_transition,
        snapshot=init_frame,
        engine=engine,
        forbidden_states=forbidden_states,
        extra_interfaces=extra_interfaces
    )
    init_conds = ratchet.run()
    init_outer = init_conds[max([s.replica for s in init_conds.samples])]
    init_inner = init_conds[0]
    return init_inner.trajectory, init_outer.trajectory


def prep_from_frame(network, initial_frame, engine):
    transitions = [trans for trans in network.sampling_transitions
                   if trans.stateA(initial_frame)]
    initial_state = transitions[0].stateA
    assert initial_state in network.all_states
    forbidden_states = set(network.all_states) - set([initial_state])
    ms_outers = network.ms_outer_objects[0]
    inners, outers = zip(*[
        bootstrap(transition, initial_frame, forbidden_states, engine,
                  ms_outers)
        for transition in transitions])

    minus_ensembles = [minus for minus in network.minus_ensembles
                       if minus.state_vol == initial_state]
    innermost = transitions[0].ensembles[0]
    trjs = [trj for trj in inners if innermost(trj)]
    assert len(trjs) > 0
    trj = trjs[0]
    minus_trajs = [
        fill_minus(minus, innermost, forbidden_states, trj, engine)
        for minus in minus_ensembles
    ]

    return minus_trajs, outers

def run_example(nsteps, setup, output_filename):
    engine, network, snapshots = setup()
    minus, outers = zip(*[prep_from_frame(network, init_frame, engine)
                          for init_frame in snapshots])
    minus = sum([list(m) for m in minus], [])
    outers = sum([list(o) for o in outers], [])

    scheme = paths.DefaultScheme(network, engine=engine)
    print("Loading initial conditions with outermost trajs")
    init_conds = scheme.initial_conditions_from_trajectories(outers)
    print("Loading initial conditions with minus trajs")
    init_conds = scheme.initial_conditions_from_trajectories(
        minus, sample_set=init_conds
    )

    storage = paths.Storage(output_filename, mode='w')
    simulation = paths.PathSampling(
        storage=storage,
        move_scheme=scheme,
        sample_set=init_conds
    )
    simulation.save_frequency = 50
    simulation.run(nsteps)
    storage.close()

def main(nsteps, examples_to_run=['mistis', 'mstis'],
         filename_format='toy_{example}_1k_OPS1_py{pyV}.nc'):
    setups = {'mstis': mstis_setup,
              'mistis': mistis_setup}
    pyV = f"{sys.version_info.major}{sys.version_info.minor}"
    for example in examples_to_run:
        print(f"Running example for {example}")
        run_example(
            nsteps=nsteps,
            setup=setups[example],
            output_filename=filename_format.format(example=example,
                                                   pyV=pyV)
        )


if __name__ == "__main__":
    # TODO: It should be possible to make a parser that handles nsteps,
    # examples to run, and filename format. Not really needed now, though.
    main(1000)

