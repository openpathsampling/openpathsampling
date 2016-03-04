from openpathsampling.features.shared import Configuration, Momentum

import mdtraj as md
import simtk.unit as u
import numpy as np
import openpathsampling as paths

import sys

__author__ = 'Jan-Hendrik Prinz'


def refresh_output(output_str, print_anyway=True, refresh=True):
    try:
        import IPython.display
    except ImportError:
        if print_anyway:
            print(output_str)
    else:
        if refresh:
            IPython.display.clear_output(wait=True)
        print(output_str)
    sys.stdout.flush()


def snapshot_from_pdb(pdb_file, simple_topology=False):
    """
    Construct a Snapshot from the first frame in a pdb file without velocities

    Parameters
    ----------
    pdb_file : str
        The filename of the .pdb file to be used

    Returns
    -------
    Snapshot
        the constructed Snapshot

    """
    pdb = md.load(pdb_file)
    velocities = np.zeros(pdb.xyz[0].shape)

    if simple_topology:
        topology = paths.Topology(*pdb.xyz[0].shape)
    else:
        topology = paths.MDTrajTopology(pdb.topology)

    configuration = Configuration(
        coordinates=u.Quantity(pdb.xyz[0], u.nanometers),
        box_vectors=u.Quantity(pdb.unitcell_vectors[0], u.nanometers),
    )

    momentum = Momentum(
        velocities=u.Quantity(velocities, u.nanometers / u.picoseconds)
    )

    snapshot = paths.Snapshot(
        topology=topology,
        configuration=configuration,
        momentum=momentum,
    )

    return snapshot


def snapshot_from_testsystem(testsystem, simple_topology=False):
    """
    Construct a Snapshot from openmm topology and state objects

    Parameters
    ----------
    omm_topology : openmm.Topology
        The filename of the .pdb file to be used

    Returns
    -------
    Snapshot
        the constructed Snapshot

    """

    velocities = u.Quantity(np.zeros(testsystem.positions.shape), u.nanometers / u.picoseconds)

    if simple_topology:
        topology = paths.Topology(*testsystem.positions.shape)
    else:
        topology = paths.MDTrajTopology(md.Topology.from_openmm(testsystem.topology))

    box_vectors = np.array([
                    v / u.nanometers for v in
                    testsystem.system.getDefaultPeriodicBoxVectors()]) * u.nanometers

    configuration = Configuration(
        coordinates=testsystem.positions,
        box_vectors=box_vectors
    )

    momentum = Momentum(
        velocities=velocities
    )

    snapshot = paths.Snapshot(
        topology=topology,
        configuration=configuration,
        momentum=momentum
    )

    return snapshot


def trajectory_from_mdtraj(mdtrajectory, simple_topology=False):
    """
    Construct a Trajectory object from an mdtraj.Trajectory object

    Parameters
    ----------
    mdtrajectory : mdtraj.Trajectory
        Input mdtraj.Trajectory

    Returns
    -------
    Trajectory
        the constructed Trajectory instance
    """

    trajectory = paths.Trajectory()
    empty_momentum = Momentum(
        velocities=u.Quantity(np.zeros(mdtrajectory.xyz[0].shape), u.nanometer / u.picosecond)
    )
    if simple_topology:
        topology = paths.Topology(*mdtrajectory.xyz[0].shape)
    else:
        topology = paths.MDTrajTopology(mdtrajectory.topology)

    for frame_num in range(len(mdtrajectory)):
        # mdtraj trajectories only have coordinates and box_vectors
        coord = u.Quantity(mdtrajectory.xyz[frame_num], u.nanometers)
        if mdtrajectory.unitcell_vectors is not None:
            box_v = u.Quantity(mdtrajectory.unitcell_vectors[frame_num],
                               u.nanometers)
        else:
            box_v = None

        config = Configuration(
            coordinates=coord,
            box_vectors=box_v
        )

        snap = paths.Snapshot(
            configuration=config,
            momentum=empty_momentum,
            topology=topology
        )
        trajectory.append(snap)

    return trajectory


def empty_snapshot_from_openmm_topology(topology, simple_topology=False):
    """
    Return an empty snapshot from an openmm.Topology object using the specified units.

    Parameters
    ----------
    topology : openmm.Topology
        the topology representing the structure and number of atoms
    units : dict of {str : simtk.unit.Unit }
        representing a dict of string representing a dimension ('length', 'velocity', 'energy') pointing the
        the simtk.unit.Unit to be used

    Returns
    -------
    Snapshot
        the complete snapshot with zero coordinates and velocities

    """
    n_atoms = topology.n_atoms

    if simple_topology:
        topology = paths.Topology(n_atoms, 3)
    else:
        topology = paths.MDTrajTopology(md.Topology.from_openmm(topology))

    configuration = Configuration(
        coordinates=u.Quantity(np.zeros((n_atoms, 3)), u.nanometers),
        box_vectors=u.Quantity(topology.setUnitCellDimensions(), u.nanometers)
    )

    momentum = Momentum(
        velocities=u.Quantity(np.zeros((n_atoms, 3)), u.nanometers / u.picoseconds)
    )

    snapshot = paths.Snapshot(
        topology=topology,
        configuration=configuration,
        momentum=momentum,
    )

    return snapshot


def to_openmm_topology(obj):
    """
    Contruct an openmm.Topology file out of a Snapshot or Configuration object. This uses the
    mdtraj.Topology in the Configuration as well as the box_vectors.

    Parameters
    ----------
    obj : Snapshot or configuration
        the object to be used in the topology construction

    Returns
    -------
    openmm.Topology
        an object representing an openmm.Topology
    """
    if obj.topology is not None:
        if hasattr(obj.topology, 'md'):
            openmm_topology = obj.topology.md.to_openmm()
            box_size_dimension = np.linalg.norm(obj.box_vectors.value_in_unit(u.nanometer), axis=1)
            openmm_topology.setUnitCellDimensions(box_size_dimension)

            return openmm_topology
    else:
        return None