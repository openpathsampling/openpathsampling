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
        this = self.__class__.__new__(self.__class__)
        AbstractSnapshot.__init__(this, topology=self.topology)
        this._is_reversed = self._is_reversed
        return this

    def create_reversed(self):
        this = self.copy()
        this._is_reversed = True
        this._reversed = self
        return this


class FeatureSnapshot(AbstractSnapshot):
    def copy(self):
        return super(FeatureSnapshot, self).copy()


@features.base.set_features(
    features.velocities,
    features.coordinates,
    features.xyz,
    features.topology
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
    features.box_vectors,
    features.xyz,
    features.topology
)
class MDSnapshot(FeatureSnapshot):
    """
    A fast MDSnapshot
    """


@features.base.set_features(
    features.configuration,
    features.momentum,
    features.xyz,
    features.topology  # for compatibility
)
class Snapshot(FeatureSnapshot):
    """
    A fast MDSnapshot
    """

