"""

@author: JD Chodera
@author: JH Prinz
"""

import abc

from openpathsampling.netcdfplus import StorableObject
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

    def __neg__(self):
        """
        Access the reversed snapshot using `-`

        Returns
        -------
        :class:`BaseSnapshot`
            the reversed copy
        """
        return self.reversed

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


def SnapshotFactory(name, features, description=None, use_lazy_reversed=False, base_class=None):
    """
    Helper to create a new Snapshot class
    
    Parameters
    ----------
    name : str
        name of the Snapshot class
    features : list of :obj:`openpathsampling.features`
        the features used to build the snapshot
    use_lazy_reversed : bool
        still in there for legacy reasons. It will make the .reversed attribute into
        a descriptor than can treat LoaderProxy objects. This feature is not relly used
        anymore and can in the best case only save little memory with slowing down construction, etc.
        Using `False` is faster
    base_class : :obj:`openpathsampling.BaseSnapshot`
        The base class the Snapshot is derived from. Default is the `BaseSnapshot` class.
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
