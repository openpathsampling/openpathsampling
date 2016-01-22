"""

@author: JD Chodera
@author: JH Prinz
"""

import mdtraj as md
import numpy as np
import simtk.unit as u
import abc

from snapshot_content import Configuration, Momentum
from openpathsampling.netcdfplus import StorableObject, lazy_loading_attributes

import features

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

    def __init__(self, topology=None, **kwargs):
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

        self._reversed = None
        self._is_reversed = False
        self.topoloty = topology

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
        if self._reversed is None:
            self._reversed = self.create_reversed()

        return self._reversed

    # ==========================================================================
    # Utility functions
    # ==========================================================================

#    @abc.abstractmethod
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
        this = self.__class__.__new__(self.__class__)
        AbstractSnapshot.__init__(this, topology=self.topology)
        this._is_reversed = self._is_reversed
        return this

    def create_reversed(self):
        this = self.copy()
        this._is_reversed = True
        this._reversed = self
        return this


# =============================================================================
# SIMULATION SNAPSHOT (COMPLETE FRAME WITH COORDINATES AND VELOCITIES)
# =============================================================================

@lazy_loading_attributes('configuration', 'momentum', '_reversed')
class Snapshot(AbstractSnapshot):
    """
    Simulation snapshot. Contains references to a configuration and momentum
    """

    __features__ = [
        features.configuration,
        features.momentum
    ]

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
        If you want to obtain a snapshot from a currently-running MD topology,
        use that topology's .current_snapshot property.

        Parameters
        ----------
        coordinates : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic coordinates (default: None)
        velocities : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic velocities (default: None)
        box_vectors : periodic box vectors (default: None)
            the periodic box vectors at current timestep (defautl: None)

        Attributes
        ----------
        coordinates : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic coordinates
        velocities : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic velocities
        box_vectors : periodic box vectors
            the periodic box vectors 
        """


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
            is_reversed=self._is_reversed,
            topology=self.topology.subset(subset)
        )
        return this


class FeatureSnapshot(AbstractSnapshot):
    def copy(self):
        return super(FeatureSnapshot, self).copy()


@features.base.set_features(
    features.velocities,
    features.coordinates,
)
class ToySnapshot(FeatureSnapshot):
    """
    Simulation snapshot. Contains references to a configuration and momentum

    Create a toy snapshot

    If you want to obtain a snapshot from a currently-running MD topology,
    use that topology's `.current_snapshot property`.
    """


@features.base.set_features(
    features.velocities,
    features.coordinates,
    features.box_vectors
)
class MDSnapshot(FeatureSnapshot):
    """
    A fast MDSnapshot
    """


@features.base.set_features(
    features.configuration,
    features.momentum,
    features.topology  # for compatibility
)
class Snapshot(FeatureSnapshot):
    """
    A fast MDSnapshot
    """

