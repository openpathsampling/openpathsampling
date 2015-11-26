__author__ = 'Jan-Hendrik Prinz'

import mdtraj as md
import simtk.unit as u
import numpy as np
import openpathsampling as paths

import sys


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


def updateunits(func):
    def inner(self, *args, **kwargs):
        my_units = {
            'length' : u.nanometer,
            'velocity' : u.nanometer / u.picoseconds,
            'energy' : u.kilojoules_per_mole
        }

        if 'units' in kwargs and kwargs['units'] is not None:
            my_units.update(kwargs['units'])

        kwargs['units'] = my_units

        return func(self, *args, **kwargs)

    return inner


@updateunits
def snapshot_from_pdb(pdb_file, units=None):
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

    snapshot = paths.Snapshot(
        coordinates=u.Quantity(pdb.xyz[0], units['length']),
        velocities=u.Quantity(velocities, units['velocity']),
        box_vectors=u.Quantity(pdb.unitcell_vectors[0], units['length']),
        potential_energy=u.Quantity(0.0, units['energy']),
        kinetic_energy=u.Quantity(0.0, units['energy']),
        topology=paths.MDTrajTopology(pdb.topology)
    )

    return snapshot

@updateunits
def snapshot_from_testsystem(testsystem, units = None):
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

    velocities = np.zeros(testsystem.positions.shape)
    topology = testsystem.topology

    box_vectors = np.array([
                    v / units['length'] for v in
                    testsystem.system.getDefaultPeriodicBoxVectors()]) * units['length']

    snapshot = paths.Snapshot(
        coordinates=testsystem.positions,
        velocities=u.Quantity(velocities, units['velocity']),
        box_vectors=box_vectors,
        potential_energy=u.Quantity(0.0, units['energy']),
        kinetic_energy=u.Quantity(0.0, units['energy']),
        topology=paths.MDTrajTopology(md.Topology.from_openmm(topology))
    )

    return snapshot

def trajectory_from_mdtraj(mdtrajectory):
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

    #TODO: Fix energies and move these to specialized CVs
    #TODO: We could also allow to have empty energies

    trajectory = paths.Trajectory()
    empty_momentum = paths.Momentum(
        velocities=u.Quantity(np.zeros(mdtrajectory.xyz[0].shape), u.nanometer / u.picosecond),
        kinetic_energy=u.Quantity(0.0, u.kilojoule_per_mole)
    )
    topology = paths.MDTrajTopology(mdtrajectory.topology)

    for frame_num in range(len(mdtrajectory)):
        # mdtraj trajectories only have coordinates and box_vectors
        coord = u.Quantity(mdtrajectory.xyz[frame_num], u.nanometers)
        if mdtrajectory.unitcell_vectors is not None:
            box_v = u.Quantity(mdtrajectory.unitcell_vectors[frame_num],
                               u.nanometers)
        else:
            box_v = None

        config = paths.Configuration(
            coordinates=coord,
            box_vectors=box_v,
            potential_energy=u.Quantity(0.0, u.kilojoule_per_mole),
            topology=topology
        )

        snap = paths.Snapshot(
            configuration=config,
            momentum=empty_momentum,
            topology=paths.MDTrajTopology(mdtrajectory.topology)
        )
        trajectory.append(snap)

    return trajectory


@updateunits
def empty_snapshot_from_openmm_topology(topology, units):
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

    snapshot = paths.Snapshot(
        coordinates=u.Quantity(np.zeros((n_atoms, 3)), units['length']),
        velocities=u.Quantity(np.zeros((n_atoms, 3)), units['velocity']),
        box_vectors=u.Quantity(topology.setUnitCellDimensions(), units['length']),
        potential_energy=u.Quantity(0.0, units['energy']),
        kinetic_energy=u.Quantity(0.0, units['energy']),
        topology=paths.MDTrajTopology(md.Topology.from_openmm(topology))
    )

    return snapshot


def units_from_snapshot(snapshot):
    """
    Returns a dict of simtk.unit.Unit instances that represent the used units in the snapshot

    Parameters
    ----------
    snapshot : Snapshot
        the snapshot to be used

    Returns
    -------
    units : dict of {str : simtk.unit.Unit }
        representing a dict of string representing a dimension ('length', 'velocity', 'energy') pointing the
        the simtk.unit.Unit to be used
    """

    units = {'velocity': None, 'energy': None, 'length': None}

    if snapshot.coordinates is not None:
        if hasattr(snapshot.coordinates, 'unit'):
            units['length'] = snapshot.coordinates.unit


    if snapshot.potential_energy is not None:
        if hasattr(snapshot.potential_energy, 'unit'):
            units['energy'] = snapshot.potential_energy.unit

    if snapshot.velocities is not None:
        if hasattr(snapshot.velocities, 'unit'):
            units['velocity'] = snapshot.velocities.unit

    return units


def to_openmm_topology(obj):
    """
    Contruct an openmm.Topology file out of a Snapshot or Configuration object. This uses the
    mdtraj.Topology in the Configuration as well as the box_vectors.

    Parameters
    ----------
    obj : Snapshot or Configuration
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
