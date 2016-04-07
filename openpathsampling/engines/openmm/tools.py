import mdtraj as md
import numpy as np
import simtk.unit as u

from snapshot import Snapshot
from topology import Topology, MDTrajTopology
from openpathsampling.engines import Trajectory

__author__ = 'Jan-Hendrik Prinz'


def snapshot_from_pdb(pdb_file, simple_topology=False):
    """
    Construct a Snapshot from the first frame in a pdb file without velocities

    Parameters
    ----------
    pdb_file : str
        The filename of the .pdb file to be used

    Returns
    -------
    :class:`openpathsampling.engines.Snapshot`
        the constructed Snapshot

    """
    pdb = md.load(pdb_file)
    velocities = np.zeros(pdb.xyz[0].shape)

    if simple_topology:
        topology = Topology(*pdb.xyz[0].shape)
    else:
        topology = MDTrajTopology(pdb.topology)

    snapshot = Snapshot.construct(
        coordinates=u.Quantity(pdb.xyz[0], u.nanometers),
        box_vectors=u.Quantity(pdb.unitcell_vectors[0], u.nanometers),
        velocities=u.Quantity(velocities, u.nanometers / u.picoseconds),
        topology=topology
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
    :class:`openpathsampling.engines.Snapshot`
        the constructed Snapshot

    """

    velocities = u.Quantity(np.zeros(testsystem.positions.shape), u.nanometers / u.picoseconds)

    if simple_topology:
        topology = Topology(*testsystem.positions.shape)
    else:
        topology = MDTrajTopology(md.Topology.from_openmm(testsystem.topology))

    box_vectors = np.array([
                    v / u.nanometers for v in
                    testsystem.system.getDefaultPeriodicBoxVectors()]) * u.nanometers

    snapshot = Snapshot.construct(
        coordinates=testsystem.positions,
        box_vectors=box_vectors,
        velocities=velocities,
        topology=topology
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
    openpathsampling.engines.Trajectory
        the constructed Trajectory instance
    """

    trajectory = Trajectory()
    empty_kinetics = Snapshot.KineticContainer(
        velocities=u.Quantity(np.zeros(mdtrajectory.xyz[0].shape), u.nanometer / u.picosecond)
    )
    if simple_topology:
        topology = Topology(*mdtrajectory.xyz[0].shape)
    else:
        topology = MDTrajTopology(mdtrajectory.topology)

    for frame_num in range(len(mdtrajectory)):
        # mdtraj trajectories only have coordinates and box_vectors
        coord = u.Quantity(mdtrajectory.xyz[frame_num], u.nanometers)
        if mdtrajectory.unitcell_vectors is not None:
            box_v = u.Quantity(mdtrajectory.unitcell_vectors[frame_num],
                               u.nanometers)
        else:
            box_v = None

        statics = Snapshot.StaticContainer(
            coordinates=coord,
            box_vectors=box_v
        )

        snap = Snapshot(
            statics=statics,
            kinetics=empty_kinetics,
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
    openpathsampling.engines.Snapshot
        the complete snapshot with zero coordinates and velocities

    """
    n_atoms = topology.n_atoms

    if simple_topology:
        topology = Topology(n_atoms, 3)
    else:
        topology = MDTrajTopology(md.Topology.from_openmm(topology))

    snapshot = Snapshot.construct(
        coordinates=u.Quantity(np.zeros((n_atoms, 3)), u.nanometers),
        box_vectors=u.Quantity(topology.setUnitCellDimensions(), u.nanometers),
        velocities=u.Quantity(np.zeros((n_atoms, 3)), u.nanometers / u.picoseconds),
        topology=topology
    )

    return snapshot


def to_openmm_topology(obj):
    """
    Contruct an openmm.Topology file out of a Snapshot or Configuration
    object. This uses the mdtraj.Topology in the Configuration as well as
    the box_vectors.

    Parameters
    ----------
    obj : openpathsampling.engines.BaseSnapshot or Configuration
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
