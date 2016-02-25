"""

@author: JD Chodera
@author: JH Prinz
"""

import abc

from openpathsampling.netcdfplus import StorableObject, lazy_loading_attributes
import features as feats


# =============================================================================
# ABSTRACT SNAPSHOT (IMPLEMENTS ONLY REVERSED SNAPSHOTS)
# =============================================================================

@lazy_loading_attributes('_reversed')
class AbstractSnapshot(StorableObject):
    """
    Simulation snapshot. Contains references to a configuration and momentum
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, topology=None):
        """
        Attributes
        ----------
        topology : openpathsamping.Topology, default: None
            The corresponding topology used with this Snapshot. Can also be None
            and means no topology is specified.
        """

        super(AbstractSnapshot, self).__init__()

        self._reversed = None
        self.topology = topology

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self is other

        return NotImplemented

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return not self is other

        return NotImplemented

    @property
    def reversed(self):
        """
        Get the reversed copy.

        Returns
        -------
        :class:`openpathsampling.snapshots.AbstractSnapshot`
            the reversed partner of the current snapshot

        Snapshots exist in pairs and this returns the reversed counter part.
        The actual implementation takes care the the reversed version have
        reversed momenta, etc. Usually these will not be stored separately but
        flipped when requested.

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

        Returns
        -------
        :class:`openpathsampling.AbstractSnapshot`
            the shallow copy

        Notes
        -----
        Shallow here means that content will not be copied but only referenced. Hence
        if you store the shallow copy it will be stored under a different idx, but the
        content (e.g. Configuration object) will not.

        """
        this = self.__class__.__new__(self.__class__)
        AbstractSnapshot.__init__(this, topology=self.topology)
        return this

    def create_reversed(self):
        this = self.copy()
        this._reversed = self
        return this


class FeatureSnapshot(AbstractSnapshot):
    def copy(self):
        return super(FeatureSnapshot, self).copy()


def snapshot_factory(name, features, description=None):
    """
    Helper to create a new Snapshot class
    
    Parameters
    ----------
    name : str
        name of the Snapshot class
    feats : 
    """

    cls = type(name, (FeatureSnapshot, ), {})
    if description is not None:
        cls.__doc__ = description

    cls = feats.set_features(*features)(cls)

    return cls


@feats.set_features(
    feats.velocities,
    feats.coordinates,
    feats.xyz,
    feats.topology
)
class ToySnapshot(FeatureSnapshot):
    """
    Simulation snapshot. Only references to coordinates and velocities
    """


@feats.set_features(
    feats.velocities,
    feats.coordinates,
    feats.box_vectors,
    feats.xyz,
    feats.topology
)
class MDSnapshot(FeatureSnapshot):
    """
    A fast MDSnapshot
    """


@feats.set_features(
    feats.configuration,
    feats.momentum,
    feats.xyz,
    feats.topology  # for compatibility
)
class Snapshot(FeatureSnapshot):
    """
    The standard MDSnapshot supporting coordinate, velocities and box_vectors
    """
