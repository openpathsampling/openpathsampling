"""

@author: JD Chodera
@author: JH Prinz
"""

import mdtraj as md
import numpy as np
import simtk.unit as u
import abc

from openpathsampling import Configuration, Momentum
from openpathsampling.netcdfplus import StorableObject, lazy_loading_attributes


def has(attr):
    def _has(func):
        def inner(self, *args, **kwargs):
            if hasattr(self, attr) and getattr(self, attr) is not None:
                return func(self, *args, **kwargs)
            else:
                return None

        return inner

    return _has


# =============================================================================
# ABSTRACT SNAPSHOT (IMPLEMENTS ONLY REVERSED SNAPSHOTS)
# =============================================================================

@lazy_loading_attributes('_reversed')
class AbstractSnapshot(StorableObject):
    """
    Simulation snapshot. Contains references to a configuration and momentum
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, is_reversed=False, reversed_copy=None, topology=None):
        """
        Attributes
        ----------
        is_reversed : bool, default: False
            The knows if in relation to potentially stored velocities these should
            be multiplied with -1 (True) or not (False).
        reversed_copy : openpathsampling.AbstractSnapshot
            If not None this is a pointer to the reversed copy (if it exists). If
            `None` (default) the reversed copy will be created and referenced.
        topology : openpathsamping.Topology, default: None
            The corresponding topology used with this Snapshot. Can also be None
            and means no topology is specified.
        """

        super(AbstractSnapshot, self).__init__()

        self.is_reversed = is_reversed
        self.topology = topology

        if reversed_copy is None:
            # this will always create the mirrored copy so we can save in pairs!
            self._reversed = self.__class__.__new__(self.__class__)
            AbstractSnapshot.__init__(
                self._reversed,
                is_reversed=not self.is_reversed,
                reversed_copy=self,
                topology=topology
            )

        else:
            self._reversed = reversed_copy

    def __eq__(self, other):
        # This implements comparison with potentially lazy loaded snapshots
        if self is other:
            return True
        elif hasattr(other, '_idx'):
            if other.__subject__ is self:
                return True

        return False

    @property
    def reversed(self):
        """
        Get the reversed copy.

        Snapshots exist in pairs and this returns the reversed counter part.
        No actual velocities are changed. Only if you ask for the velocities of
        a reversed object the velocities will be multiplied by -1.
        """
        return self._reversed

    # ==========================================================================
    # Utility functions
    # ==========================================================================

    @abc.abstractmethod
    def copy(self):
        """
        Returns a shallow copy of the instance itself. The contained
        configuration and momenta are not copied.

        This will also create a new reversed copy automatically.

        Returns
        -------
        openpathsampling.AbstractSnapshot()
            the shallow object

        Notes
        -----
        Shallow here means that content will not be copied but only referenced. Hence
        if you store the shallow copy it will be stored under a different idx, but the
        content (e.g. Configuration object) will not.

        """
        this = AbstractSnapshot(is_reversed=self.is_reversed, topology=self.topology)
        return this

    def reversed_copy(self):
        """
        Returns a shallow reversed copy of the instance itself. The
        contained configuration and momenta are not copied and the momenta
        are marked reversed.

        Since Snapshots are made in pairs this will also lead to a new
        (non-)reversed copy!

        Returns
        -------
        Snapshot()
            the reversed shallow copy
        """

        obj = self.copy()
        return obj.reversed


# =============================================================================
# SIMULATION SNAPSHOT (COMPLETE FRAME WITH COORDINATES AND VELOCITIES)
# =============================================================================

@lazy_loading_attributes('configuration', 'momentum', '_reversed')
class Snapshot(AbstractSnapshot):
    """
    Simulation snapshot. Contains references to a configuration and momentum
    """

    __features__ = ['Configurations', 'Momenta']

    def __init__(self, coordinates=None, velocities=None, box_vectors=None,
                 potential_energy=None, kinetic_energy=None, topology=None,
                 configuration=None, momentum=None, is_reversed=False,
                 reversed_copy=None):
        """
        Create a simulation snapshot. Initialization happens primarily in
        one of two ways:
            1. Specify `Configuration` and `Momentum` objects
            2. Specify the things which make up `Configuration` and
               `Momentum` objects, i.e., coordinates, velocities, box
               vectors, etc.
        If you want to obtain a snapshot from a currently-running MD engine,
        use that engine's .current_snapshot property.

        Parameters
        ----------
        coordinates : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic coordinates (default: None)
        velocities : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic velocities (default: None)
        box_vectors : periodic box vectors (default: None)
            the periodic box vectors at current timestep (defautl: None)
        potential_energy : simtk.unit.Quantity of units energy/mole
            potential energy at current timestep (default: None)
        kinetic_energy : simtk.unit.Quantity of units energy/mole
            kinetic energy at current timestep (default: None)
            
        Attributes
        ----------
        coordinates : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic coordinates
        velocities : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic velocities
        box_vectors : periodic box vectors
            the periodic box vectors 
        potential_energy : simtk.unit.Quantity of units energy/mole
            potential energy
        kinetic_energy : simtk.unit.Quantity of units energy/mole
            kinetic energy
        """

        super(Snapshot, self).__init__(is_reversed, reversed_copy, topology)

        if configuration is None and momentum is None:
            if coordinates is not None:
                configuration = Configuration(
                    coordinates=coordinates,
                    box_vectors=box_vectors,
                    potential_energy=potential_energy
                )

            if velocities is not None:
                momentum = Momentum(
                    velocities=velocities,
                    kinetic_energy=kinetic_energy
                )

        self.configuration = configuration
        self.momentum = momentum

        if reversed_copy is None:
            self._reversed.configuration = self.configuration
            self._reversed.momentum = self.momentum

    @property
    @has('configuration')
    def coordinates(self):
        """
        The coordinates in the configuration
        """
        return self.configuration.coordinates

    @property
    @has('momentum')
    def velocities(self):
        """
        The velocities in the configuration. If the snapshot is reversed a
        copy of the original (unreversed) velocities is made which is then
        returned
        """
        if self.is_reversed:
            return -1.0 * self.momentum.velocities
        else:
            return self.momentum.velocities

    @property
    @has('configuration')
    def box_vectors(self):
        """
        The box_vectors in the configuration
        """
        return self.configuration.box_vectors

    @property
    @has('configuration')
    def potential_energy(self):
        """
        The potential_energy in the configuration
        """
        return self.configuration.potential_energy

    @property
    @has('momentum')
    def kinetic_energy(self):
        """
        The kinetic_energy in the momentum
        """
        return self.momentum.kinetic_energy

    @property
    @has('configuration')
    def n_atoms(self):
        """
        The number of atoms in the snapshot
        """
        return self.coordinates.shape[0]

    @property
    @has('configuration')
    @has('momentum')
    def total_energy(self):
        """
        The total energy (sum of potential and kinetic) of the snapshot
        """
        return self.kinetic_energy + self.potential_energy

    # ==========================================================================
    # Utility functions
    # ==========================================================================

    @property
    @has('configuration')
    def xyz(self):
        """
        Coordinates without dimensions.

        SERIOUS PROBLEM: whether .coordinates returns a u.Quantity or jut
        numpy array depending on situation (at least for ToyDynamics). This
        is bad.
        """
        coord = self.configuration.coordinates
        if type(coord) is u.Quantity:
            return coord._value
        else:
            return coord

    def copy(self):
        this = self.__class__(
            configuration=self.configuration,
            momentum=self.momentum,
            is_reversed=self.is_reversed,
            topology=self.topology
        )
        return this

    @has('configuration')
    def md(self):
        """
        Returns a mdtraj Trajectory object that contains only one frame

        Returns
        -------
        mdtraj.Trajectory
            the actual trajectory object. Can be used with all functions from mdtraj

        Notes
        -----
        Rather slow since the topology has to be made each time. Try to avoid it
        """

        n_atoms = self.n_atoms

        output = np.zeros([1, n_atoms, 3], np.float32)
        output[0, :, :] = self.coordinates

        return md.Trajectory(output, self.topology.md)

    def subset(self, subset):
        """
        Return a deep copy of the snapshot with reduced set of coordinates. Takes also care
        of adjusting the topology.

        Parameters
        ----------
        subset : list of int
            a list of atomic indices specifying which entries to keep.

        Notes
        -----
        So far the potential and kinetic energies are copied and are thus false but still useful!?!

        """

        this = Snapshot(
            configuration=self.configuration.copy(subset),
            momentum=self.momentum.copy(subset),
            is_reversed=self.is_reversed,
            topology=self.topology.subset(subset)
        )
        return this


@lazy_loading_attributes('_reversed')
class ToySnapshot(AbstractSnapshot):
    """
    Simulation snapshot. Contains references to a configuration and momentum
    """

    __features__ = ['Velocities', 'Coordinates']

    def __init__(self, coordinates=None, velocities=None, is_reversed=False, topology=None,
                 reversed_copy=None):
        """
        Create a toy snapshot

        If you want to obtain a snapshot from a currently-running MD engine,
        use that engine's `.current_snapshot property`.

        Parameters
        ----------
        coordinates : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic coordinates (default: None)
        velocities : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic velocities (default: None)

        Attributes
        ----------
        coordinates : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic coordinates
        velocities : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic velocities
        """

        super(ToySnapshot, self).__init__(is_reversed, reversed_copy, topology)

        self.coordinates = coordinates
        self.velocities = velocities

        if reversed_copy is None:
            self._reversed.coordinates = self.coordinates
            self._reversed.velocities = -1.0 * self.velocities

    # ==========================================================================
    # Utility functions
    # ==========================================================================

    @property
    def xyz(self):
        """
        Coordinates without dimensions.

        """
        return self.coordinates

    def copy(self):
        this = ToySnapshot(
            self.coordinates,
            self.velocities,
            is_reversed=self.is_reversed,
            topology=self.topology
        )
        return this
