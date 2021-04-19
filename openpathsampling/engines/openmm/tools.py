import numpy as np

from openpathsampling.integration_tools import (
    md, error_if_no_mdtraj, unit, error_if_no_simtk_unit
)
try:
    import simtk.openmm
    import simtk.openmm.app
    from simtk import unit as u
    from simtk.openmm.app.internal.unitcell import reducePeriodicBoxVectors
except ImportError:
    # this happens when we directly import tools (e.g., for
    # trajectory_to/from_mdtraj) when we don't have OpenMM installed. In
    # that case, we skip the imports in the else statement
    # (engines/openmm/__init__.py prevents them from being made)
    pass
else:
    from .snapshot import Snapshot
    from openpathsampling.engines.topology import Topology, MDTrajTopology

from openpathsampling.engines import Trajectory, NoEngine, SnapshotDescriptor


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
    error_if_no_simtk_unit("snapshot_from_testsystem")
    u_nm = unit.nanometers
    u_ps = unit.picoseconds
    velocities = unit.Quantity(np.zeros(testsystem.positions.shape),
                               u_nm / u_ps)

    if simple_topology:
        topology = Topology(*testsystem.positions.shape)
    else:
        topology = MDTrajTopology(md.Topology.from_openmm(testsystem.topology))

    if periodic:
        sys_box_vectors = testsystem.system.getDefaultPeriodicBoxVectors()
        box_vectors = np.array([v / u_nm for v in sys_box_vectors]) * u_nm
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
    error_if_no_simtk_unit("trajectory_from_mdtraj")
    trajectory = Trajectory()
    u_nm = unit.nanometer
    u_ps = unit.picosecond
    vel_unit = u_nm / u_ps

    if simple_topology:
        topology = Topology(*mdtrajectory.xyz[0].shape)
    else:
        topology = MDTrajTopology(mdtrajectory.topology)

    if velocities is None:
        empty_vel = unit.Quantity(np.zeros(mdtrajectory.xyz[0].shape),
                                  vel_unit)

    if mdtrajectory.unitcell_vectors is not None:
        box_vects = unit.Quantity(mdtrajectory.unitcell_vectors,
                                  unit.nanometers)
    else:
        box_vects = [None] * len(mdtrajectory)


    engine = TopologyEngine(topology)

    for frame_num in range(len(mdtrajectory)):
        # mdtraj trajectories only have coordinates and box_vectors
        coord = unit.Quantity(mdtrajectory.xyz[frame_num], u_nm)
        if velocities is not None:
            vel = unit.Quantity(velocities[frame_num], vel_unit)
        else:
            vel = empty_vel

        box_v = box_vects[frame_num]

        statics = Snapshot.StaticContainer(
            coordinates=coord,
            box_vectors=box_v,
            engine=engine
        )
        kinetics = Snapshot.KineticContainer(velocities=vel,
                                             engine=engine)

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

    error_if_no_simtk_unit("empty_snapshot_from_openmm_topology")
    u_nm = unit.nanometers
    u_ps = unit.picoseconds
    n_atoms = topology.n_atoms

    if simple_topology:
        topology = Topology(n_atoms, 3)
    else:
        error_if_no_mdtraj("empty_snaphsot_from_openmm_topology")
        topology = MDTrajTopology(md.Topology.from_openmm(topology))

    snapshot = Snapshot.construct(
        coordinates=unit.Quantity(np.zeros((n_atoms, 3)), u_nm),
        box_vectors=unit.Quantity(topology.setUnitCellDimensions(), u_nm),
        velocities=unit.Quantity(np.zeros((n_atoms, 3)), u_nm / u_ps),
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
                obj.box_vectors.value_in_unit(unit.nanometer), axis=1)
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

    # For now, let's keep all the code in one place, and better for
    # engines.openmm.tools to require engines.trajectory than vice versa
    return trajectory.to_mdtraj(md_topology)

def ops_load_trajectory(filename, **kwargs):
    error_if_no_mdtraj("ops_load_trajectory")
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
    nm = unit.nanometer
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

def n_dofs_from_system(system):
    """Get the number of degrees of freedom from an Openmm System

    Parameters
    ----------
    system : :class:`simtk.openmm.System`
        object describing the system

    Returns
    -------
    int :
        number of degrees of freedom
    """
    # dof calculation based on OpenMM's StateDataReporter
    n_spatial = 3
    n_particles = system.getNumParticles()
    dofs_particles = sum([n_spatial for i in range(n_particles)
                          if system.getParticleMass(i) > 0*u.dalton])
    dofs_constaints = system.getNumConstraints()
    dofs_motion_removers = 0
    has_cm_motion_remover = any(
        type(system.getForce(i)) == simtk.openmm.CMMotionRemover
        for i in range(system.getNumForces())
    )
    if has_cm_motion_remover:
        dofs_motion_removers += 3
    dofs = dofs_particles - dofs_constaints - dofs_motion_removers
    return dofs

def has_constraints_from_system(system):
    """Get the number of degrees of freedom from an Openmm System

    Parameters
    ----------
    system : :class:`simtk.openmm.System`
        object describing the system

    Returns
    -------
    bool :
        whether there are constraints
    """
    return system.getNumConstraints() > 0
