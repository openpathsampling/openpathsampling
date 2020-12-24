"""

@author: JD Chodera
@author: JH Prinz
"""

import abc

from openpathsampling.netcdfplus import StorableObject
from . import features as feats


# =============================================================================
# ABSTRACT SNAPSHOT (IMPLEMENTS ONLY REVERSED SNAPSHOTS)
# =============================================================================

class BaseSnapshot(StorableObject):
    """
    Simulation snapshot. Contains references to a configuration and momentum

    Parameters
    ----------
    topology : openpathsamping.Topology, default: None
        The corresponding topology used with this Snapshot. Can also be None
        and means no topology is specified.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, topology=None):
        super(BaseSnapshot, self).__init__()

        self._reversed = None
        self.topology = topology

    # def __eq__(self, other):
    #     if self is other:
    #         return True
    #
    #     if isinstance(other, BaseSnapshot):
    #         return self.__uuid__ == other.__uuid__
    #
    #     return NotImplemented
    #
    # def __ne__(self, other):
    #     return not self == other

    # def __hash__(self):
    #     return self.__uuid__ & 1152921504606846975

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
        Shallow here means that content will not be copied but only referenced.
        Hence if you store the shallow copy it will be stored under a different
        idx, but the content (e.g. Configuration object) will not.

        """
        this = self.__class__.__new__(self.__class__)
        BaseSnapshot.__init__(this, topology=self.topology)
        return this

    def copy_with_replacement(self, **kwargs):
        cp = self.copy()  # this will copy all, but it is simple
        for key, value in kwargs.items():
            if hasattr(cp, key):
                setattr(cp, key, value)
            else:
                raise TypeError("copy_with_replacement() got an "
                                "unexpected keyword argument '%s'" % key)

        return cp

    def create_reversed(self):
        this = self.copy()
        this._reversed = self
        return this


def SnapshotFactory(
        name,
        features,
        description=None,
        use_lazy_reversed=False,
        base_class=None):
    """
    Helper to create a new Snapshot class

    Parameters
    ----------
    name : str
        name of the Snapshot class
    features : list of :obj:`openpathsampling.features`
        the features used to build the snapshot
    description : str
        the string to be used as basis for the docstring of the new class
        it will be merged with the docs for the features
    use_lazy_reversed : bool
        still in there for legacy reasons. It will make the .reversed attribute
        into a descriptor than can treat LoaderProxy objects. This feature is
        not relly used anymore and can in the best case only save little memory
        with slowing down construction, etc. Using `False` is faster
    base_class : :obj:`openpathsampling.BaseSnapshot`
        The base class the Snapshot is derived from.
        Default is the `BaseSnapshot` class.

    Returns
    -------
    :class:`openpathsampling.Snapshot`
        the created `Snapshot` class
    """

    if base_class is None:
        base_class = BaseSnapshot

    if type(base_class) is not tuple:
        base_class = (base_class,)

    cls = type(name, base_class, {})
    if description is not None:
        cls.__doc__ = description

    cls = feats.attach_features(
        features,
        use_lazy_reversed=use_lazy_reversed)(cls)

    return cls


class SnapshotDescriptor(frozenset, StorableObject):
    """Container for information about snapshots generated by an engine.

    Snapshot descriptors are used to define the dimensions of the features
    used in a snapshot, in order to set correct sizes in storage. For
    example, the arrays of atomic positions and velocities will each be of
    shape ``(n_atoms, n_spatial)``. The snapshot descriptor stores the
    values of ``n_atoms`` and ``n_spatial``. It also knows the class of
    snapshot to be created by the engine. This is usually created upon
    initialization of the engine, using information from the engine's
    initialization parameters.

    In practice, it is probably easiest to create snapshot descriptors using
    their :meth:`.construct` method.

    Parameters
    ----------
    contents : list of 2-tuples
        Key-value pairs for information to be stored. One must be the key
        'class', mapped to a snapshot class.
    """
    def __init__(self, contents):
        StorableObject.__init__(self)
        frozenset.__init__(contents)
        self._dimensions = dict(self)
        self._cls = self._dimensions['class']
        del self._dimensions['class']

    @property
    def snapshot_class(self):
        return self._cls

    @property
    def dimensions(self):
        return self._dimensions

    @classmethod
    def from_dict(cls, dct):
        return cls(dct.items())

    def to_dict(self):
        return dict(self)

    @staticmethod
    def construct(snapshot_class, snapshot_dimensions):
        """Convenience method to create a snapshot descriptor.

        Parameters
        ----------
        snapshot_class : class
            class that creates snapshots for this engine
        snapshot_dimensions : dict
            dictionary mapping dimension name to integer

        Returns
        -------
        :class:`.SnapshotDescriptor`:
            descriptor based on the input information

        Examples
        --------
        >>> from openpathsampling.engines import SnapshotDescriptor, toy
        >>> descriptor = SnapshotDescriptor.construct(
        ...     snapshot_class=toy.Snapshot,
        ...     snapshot_dimensions={'n_atoms': 1, 'n_spatial': 2}
        ... )

        """
        d = {'class': snapshot_class}

        if set(snapshot_class.__features__.dimensions) > \
                set(snapshot_dimensions.keys()):
            raise RuntimeError(
                ('Snapshot of type %s needs %s as dimensions, '
                 'you only provided %s') %
                (
                    snapshot_class.__class__.__name__,
                    str(set(snapshot_class.__features__['dimensions'])),
                    str(set(snapshot_dimensions.keys()))
                )
            )

        d.update(snapshot_dimensions)

        return SnapshotDescriptor(list(d.items()))
