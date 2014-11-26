__author__ = 'jan-hendrikprinz'

import mdtraj as md
from opentis.snapshot import Configuration, Snapshot, Momentum
from opentis.trajectory import Trajectory
import simtk.unit as u
import numpy as np




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
def configuration_from_pdb(pdb_file, units = None):
    pdb = md.load(pdb_file)

    configuration = Configuration(
        coordinates=u.Quantity(pdb.xyz[0], units['length']),
        box_vectors=u.Quantity(pdb.unitcell_vectors, units['length']),
        potential_energy=u.Quantity(0.0, units['energy']),
        topology=pdb.topology
    )

    return configuration

@updateunits
def snapshot_from_pdb(pdb_file, units = None):
    pdb = md.load(pdb_file)

    velocities = np.zeros(pdb.xyz[0].shape)

    snapshot = Snapshot(
        coordinates=u.Quantity(pdb.xyz[0], units['length']),
        velocities=u.Quantity(velocities, units['velocity']),
        box_vectors=u.Quantity(pdb.unitcell_vectors, units['length']),
        potential_energy=u.Quantity(0.0, units['energy']),
        kinetic_energy=u.Quantity(0.0, units['energy']),
        topology=pdb.topology
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
    """
    trajectory = Trajectory()
    empty_momentum = Momentum()
    for frame_num in range(mdtrajectory.n_frames):
        # mdtraj trajectories only have coordinates and box_vectors
        coord = u.Quantity(mdtrajectory.xyz[frame_num], u.nanometers)
        if mdtrajectory.unitcell_vectors is not None:
            box_v = u.Quantity(mdtrajectory.unitcell_vectors[frame_num],
                             u.nanometers)
        else:
            box_v = None
        config = Configuration(coordinates=coord, box_vectors=box_v)

        snap = Snapshot(configuration=config, momentum=empty_momentum, topology=mdtrajectory.topology)
        trajectory.append(snap)

    return trajectory

@updateunits
def empty_snapshot_from_openmm_topology(topology, units):
    n_atoms = topology.n_atoms

    snapshot = Snapshot(
        coordinates=u.Quantity(np.zeros((n_atoms, 3)), units['length']),
        velocities=u.Quantity(np.zeros((n_atoms, 3)), units['velocity']),
        box_vectors=u.Quantity(topology.setUnitCellDimensions(), units['length']),
        potential_energy=u.Quantity(0.0, units['energy']),
        kinetic_energy=u.Quantity(0.0, units['energy']),
        topology=md.Topology.from_openmm(topology)
    )

    return snapshot

def units_from_snapshot(snapshot):
    return {
        'length' : snapshot.coordinates.unit,
        'velocity' : snapshot.velocities.unit,
        'energy' : snapshot.potential_energy.unit
    }

def to_openmm_topology(obj):
    if obj.topology is not None:
        openmm_topology = obj.topology.to_openmm()
        box_size_dimension = np.linalg.norm(obj.box_vectors.value_in_unit(u.nanometer), axis=1)[0]
        openmm_topology.setUnitCellDimensions(box_size_dimension)

        return openmm_topology
    else:
        return None
