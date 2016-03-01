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

class BaseSnapshot(StorableObject):
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

        super(BaseSnapshot, self).__init__()

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

    def copy(self):
        """
        Returns a shallow copy of the instance itself. The contained
        configuration and momenta are not copied.

        Returns
        -------
        :class:`openpathsampling.BaseSnapshot`
            the shallow copy

        Notes
        -----
        Shallow here means that content will not be copied but only referenced. Hence
        if you store the shallow copy it will be stored under a different idx, but the
        content (e.g. Configuration object) will not.

        """
        this = self.__class__.__new__(self.__class__)
        BaseSnapshot.__init__(this, topology=self.topology)
        return this

    def create_reversed(self):
        this = self.copy()
        this._reversed = self
        return this


def SnapshotFactory(name, features, description=None, use_lazy_reversed=True, base_class=None):
    """
    Helper to create a new Snapshot class
    
    Parameters
    ----------
    name : str
        name of the Snapshot class
    features : list of :obj:`openpathsampling.features`
        the features used to build the snapshot
    use_lazy_reversed : bool
    base_class : :obj:`openpathsampling.AbstractSnapshot`
    """


    if base_class is None:
        base_class = BaseSnapshot

    if type(base_class) is not tuple:
        base_class = (base_class, )

    cls = type(name, base_class, {})
    if description is not None:
        cls.__doc__ = description

    cls = feats.attach_features(features, use_lazy_reversed=use_lazy_reversed)(cls)

    return cls


@feats.attach_features([
    feats.velocities,
    feats.coordinates,
    feats.xyz,
    feats.topology
])
class ToySnapshot(BaseSnapshot):
    """
    Simulation snapshot. Only references to coordinates and velocities
    """


@feats.attach_features([
    feats.velocities,
    feats.coordinates,
    feats.box_vectors,
    feats.xyz,
    feats.topology
])
class MDSnapshot(BaseSnapshot):
    """
    A fast MDSnapshot
    """


@feats.attach_features([
    feats.configuration,
    feats.momentum,
    feats.xyz,
    feats.topology  # for compatibility
])
class Snapshot(BaseSnapshot):
    """
    The standard MDSnapshot supporting coordinate, velocities and box_vectors
    """
