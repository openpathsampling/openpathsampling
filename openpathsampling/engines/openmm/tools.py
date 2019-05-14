import mdtraj as md
import numpy as np
import simtk.unit as u

from .snapshot import Snapshot
from .topology import Topology, MDTrajTopology
from openpathsampling.engines import Trajectory, NoEngine, SnapshotDescriptor

from simtk.openmm.app.internal.unitcell import reducePeriodicBoxVectors

__author__ = 'Jan-Hendrik Prinz'


class TopologyEngine(NoEngine):
    _default_options = {}

    def __init__(self, topology):

        descriptor = SnapshotDescriptor.construct(
            Snapshot,
            {
                'n_atoms': topology.n_atoms,
                'n_spatial': topology.n_spatial
            }
        )

        super(NoEngine, self).__init__(
            descriptor=descriptor
        )

        self.topology = topology

    @property
    def mdtraj_topology(self):
        return self.topology.mdtraj

    def to_dict(self):
        return {
            'topology': self.topology,
        }


class FileEngine(TopologyEngine):

    _default_options = {}

    def __init__(self, topology, filename):
        super(FileEngine, self).__init__(
            topology=topology
        )

        self.filename = filename

    def to_dict(self):
        return {
            'topology': self.topology,
            'filename': self.filename
        }


class OpenMMToolsTestsystemEngine(TopologyEngine):

    _default_options = {}

    def __init__(self, topology, testsystem_name):

        super(OpenMMToolsTestsystemEngine, self).__init__(
            topology=topology
        )

        self.testsystem_name = testsystem_name

    def to_dict(self):
        return {
            'topology': self.topology,
            'testsystem_name': self.testsystem_name
        }


def snapshot_from_pdb(pdb_file, simple_topology=False):
    """
    Construct a Snapshot from the first frame in a pdb file without velocities

    Parameters
    ----------
    pdb_file : str
        The filename of the .pdb file to be used
    simple_topology : bool
        if `True` only a simple topology with n_atoms will be created.
        This cannot be used with complex CVs but loads and stores very fast

    Returns
    -------
    :class:`openpathsampling.engines.Snapshot`
        the constructed Snapshot

    """
    snap = ops_load_trajectory(pdb_file)[0]

    if simple_topology:
        topology = Topology(*pdb.xyz[0].shape)
    else:
        topology = snap.topology

    snapshot = Snapshot.construct(
        coordinates=snap.coordinates,
        box_vectors=snap.box_vectors,
        velocities=snap.velocities,
        engine=FileEngine(topology, pdb_file)
    )

    return snapshot


def topology_from_pdb(pdb_file, simple_topology=False):
    """
    Construct a Topology from the first frame in a pdb file without velocities

    Parameters
    ----------
    pdb_file : str
        The filename of the .pdb file to be used
    simple_topology : bool
        if `True` only a simple topology with n_atoms will be created.
        This cannot be used with complex CVs but loads and stores very fast

    Returns
    -------
    :class:`openpathsampling.engines.Snapshot`
        the constructed Snapshot

    """
    pdb = md.load(pdb_file)

    if simple_topology:
        topology = Topology(*pdb.xyz[0].shape)
    else:
        topology = MDTrajTopology(pdb.topology)

    return topology


def snapshot_from_testsystem(testsystem, simple_topology=False,
                             periodic=True):
    """
    Construct a Snapshot from openmm topology and state objects

    Parameters
    ----------
    testsystem : openmmtools.Topology
        The filename of the .pdb file to be used
    simple_topology : bool
        if `True` only a simple topology with n_atoms will be created.
        This cannot be used with complex CVs but loads and stores very fast
    periodic : bool
        True (default) if system is periodic; if False, box vectors are None

    Returns
    -------
    :class:`openpathsampling.engines.Snapshot`
        the constructed Snapshot

    """

    velocities = u.Quantity(
        np.zeros(testsystem.positions.shape), u.nanometers / u.picoseconds)

    if simple_topology:
        topology = Topology(*testsystem.positions.shape)
    else:
        topology = MDTrajTopology(md.Topology.from_openmm(testsystem.topology))

    if periodic:
        box_vectors = \
            np.array([
                v / u.nanometers for v in
                testsystem.system.getDefaultPeriodicBoxVectors()]) * u.nanometers
    else:
        box_vectors = None

    snapshot = Snapshot.construct(
        coordinates=testsystem.positions,
        box_vectors=box_vectors,
        velocities=velocities,
        engine=OpenMMToolsTestsystemEngine(topology, testsystem.name)
    )

    return snapshot


def trajectory_from_mdtraj(mdtrajectory, simple_topology=False,
                           velocities=None):
    """
    Construct a Trajectory object from an mdtraj.Trajectory object

    Parameters
    ----------
    mdtrajectory : mdtraj.Trajectory
        Input mdtraj.Trajectory
    simple_topology : bool
        if `True` only a simple topology with n_atoms will be created.
        This cannot be used with complex CVs but loads and stores very fast
    velocities : np.array
        velocities in units of nm/ps

    Returns
    -------
    openpathsampling.engines.Trajectory
        the constructed Trajectory instance
    """
    trajectory = Trajectory()
    vel_unit = u.nanometer / u.picosecond

    if simple_topology:
        topology = Topology(*mdtrajectory.xyz[0].shape)
    else:
        topology = MDTrajTopology(mdtrajectory.topology)

    if velocities is None:
        empty_vel = u.Quantity(np.zeros(mdtrajectory.xyz[0].shape),
                               vel_unit)


    engine = TopologyEngine(topology)

    for frame_num in range(len(mdtrajectory)):
        # mdtraj trajectories only have coordinates and box_vectors
        coord = u.Quantity(mdtrajectory.xyz[frame_num], u.nanometers)
        if velocities is not None:
            vel = u.Quantity(velocities[frame_num], vel_unit)
        else:
            vel = empty_vel

        if mdtrajectory.unitcell_vectors is not None:
            box_v = u.Quantity(mdtrajectory.unitcell_vectors[frame_num],
                               u.nanometers)
        else:
            box_v = None

        statics = Snapshot.StaticContainer(
            coordinates=coord,
            box_vectors=box_v
        )
        kinetics = Snapshot.KineticContainer(velocities=vel)

        snap = Snapshot(
            statics=statics,
            kinetics=kinetics,
            engine=engine
        )
        trajectory.append(snap)

    return trajectory


def empty_snapshot_from_openmm_topology(topology, simple_topology=False):
    """
    Return an empty snapshot from an openmm.Topology object

    Velocities will be set to zero.

    Parameters
    ----------
    topology : openmm.Topology
        the topology representing the structure and number of atoms
    simple_topology : bool
        if `True` only a simple topology with n_atoms will be created.
        This cannot be used with complex CVs but loads and stores very fast

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
        velocities=u.Quantity(
            np.zeros((n_atoms, 3)), u.nanometers / u.picoseconds),
        engine=TopologyEngine(topology)
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
        if hasattr(obj.topology, 'mdtraj'):
            openmm_topology = obj.topology.mdtraj.to_openmm()
            box_size_dimension = np.linalg.norm(
                obj.box_vectors.value_in_unit(u.nanometer), axis=1)
            openmm_topology.setUnitCellDimensions(box_size_dimension)

            return openmm_topology
    else:
        return None


def trajectory_to_mdtraj(trajectory, md_topology=None):
    """
    Construct a `mdtraj.Trajectory` object from an :obj:`Trajectory` object

    Parameters
    ----------
    trajectory : :obj:`openpathsampling.engines.Trajectory`
        Input Trajectory

    Returns
    -------
    :obj:`mdtraj.Trajectory`
        the constructed Trajectory instance
    """
    if not hasattr(trajectory, 'to_mdtraj'):
        try:
            _ = len(trajectory)
        except TypeError:
            trajectory = Trajectory([trajectory])
        else:
            trajectory = Trajectory(trajectory)

    # TODO: The following would work if we remove trajectory.to_mdtraj()
    # For now, let's keep all the code in one place, and better for
    # engines.openmm.tools to require engines.trajectory than vice versa
    # output = trajectory.xyz
    # traj = md.Trajectory(output, md_topology)
    # traj.unitcell_vectors = trajectory.box_vectors
    return trajectory.to_mdtraj(md_topology)

def ops_load_trajectory(filename, **kwargs):
    return trajectory_from_mdtraj(md.load(filename, **kwargs))

# ops_load_trajectory and the mdtraj stuff is not OpenMM-specific

def reduced_box_vectors(snapshot):
    """Reduced box vectors for a snapshot (with units)

    See also
    --------
    reduce_trajectory_box_vectors

    Parameters
    ----------
    snapshot : :class:`.Snapshot`
        input snapshot

    Returns
    -------
    :class:`.Snapshot`
        snapshot with correctly reduced box vectors
    """
    nm = u.nanometer
    return np.array(
        reducePeriodicBoxVectors(snapshot.box_vectors).value_in_unit(nm)
    ) * nm

def reduce_trajectory_box_vectors(trajectory):
    """Trajectory with reduced box vectors.

    OpenMM has strict requirements on the box vectors describing the unit
    cell. In some cases, such as trajectories loaded from files that have
    rounded the box vectors, these conditions might not be satisfied. This
    method forces the box vectors to meet OpenMM's criteria.

    Parameters
    ----------
    trajectory : :class:`.Trajectory`
        input trajectory

    Returns
    -------
    :class:`.Trajectory`
        trajectory with correctly reduced box vectors
    """
    return Trajectory([
        snap.copy_with_replacement(box_vectors=reduced_box_vectors(snap))
        for snap in trajectory
    ])


def load_trr(trr_file, top, velocities=True):
    """Load a TRR file, ready for use as input to an OpenMMEngine.

    This is a single method to handle several peculiarities of both the TRR
    format (which rounds some values) and OpenMM (which has certain
    requirements of box vectors), plus the possibility that you'll want
    velocities.

    Parameters
    ----------
    trr_file : string
        name of TRR file to load
    top : string
        name of topology (e.g., ``.gro``) file to use. See MDTraj
        documentation on md.load.
    velocities : bool
        whether to also load velocities from the TRR file; default ``True``

    Return
    ------
    :class:`.Trajectory`
        the OPS trajectory, with OpenMM-reduced box vectors and velocities
        (if requested)
    """
    mdt = md.load(trr_file, top=top)
    trr = md.formats.TRRTrajectoryFile(trr_file)
    if velocities:
        vel = trr._read(n_frames=len(mdt), atom_indices=None,
                        get_velocities=True)[5]
    else:
        vel = None

    traj = trajectory_from_mdtraj(mdt, velocities=vel)
    return reduce_trajectory_box_vectors(traj)
