from openpathsampling.snapshot import Snapshot, AbstractSnapshot, ToySnapshot
from openpathsampling.trajectory import Trajectory
from openpathsampling.netcdfplus import ObjectStore, LoaderProxy

import snapshot_features as ft
from snapshot_features import ConfigurationStore, MomentumStore


# =============================================================================================
# ABSTRACT BASE CLASS FOR SNAPSHOTS
# =============================================================================================

class AbstractSnapshotStore(ObjectStore):
    """
    An ObjectStore for Snapshots in netCDF files.
    """

    def __init__(self, snapshot_class):
        """

        Attributes
        ----------
        snapshot_class : openpathsampling.AbstractSnapshot
            a snapshot class that this Store is supposed to store

        """
        super(AbstractSnapshotStore, self).__init__(AbstractSnapshot, json=False)
        self.snapshot_class = snapshot_class

    def to_dict(self):
        return {
            'snapshot_class': self.snapshot_class
        }

    def _get(self, idx, from_reversed=False):
        if from_reversed:
            obj = self.cache[AbstractSnapshotStore.paired_idx(idx)]

            return AbstractSnapshot(
                is_reversed=not obj.is_reversed,
                topology=obj.topology,
                reversed_copy=LoaderProxy(self, AbstractSnapshotStore.paired_idx(idx))
            )
        else:
            momentum_reversed = self.vars['momentum_reversed'][idx]
            topology = self.storage.topology

            return AbstractSnapshot(
                is_reversed=momentum_reversed,
                topology=topology,
                reversed_copy=LoaderProxy(self, AbstractSnapshotStore.paired_idx(idx))
            )

    def _load(self, idx):
        """
        Load a snapshot from the storage.

        Parameters
        ----------
        idx : int
            the integer index of the snapshot to be loaded

        Returns
        -------
        snapshot : Snapshot
            the loaded snapshot instance
        """

        try:
            return self._get(idx, True)
        except KeyError:
            return self._get(idx)

    @staticmethod
    def paired_idx(idx):
        """
        Return the paired index

        Snapshots are stored in pairs (2n, 2n+1) where one is the reversed copy.
        This make storing CVs easier. This function allows to get the paired index
        or the index of snapshot.reversed

        The implementation uses the trick that all you have to do is flip the lowest bit
        that determines even or odd.

        Parameters
        ----------
        idx : int
            the one part of the paired index

        Returns
        -------
        int
            the other part of the paired index
        """
        return idx ^ 1

    def _set(self, idx, snapshot):
        self.vars['momentum_reversed'][idx] = snapshot.is_reversed
        self.vars['momentum_reversed'][AbstractSnapshotStore.paired_idx(idx)] = not snapshot.is_reversed

    def _save(self, snapshot, idx):
        """
        Add the current state of the snapshot in the database.

        Parameters
        ----------
        snapshot : Snapshot()
            the snapshot to be saved
        idx : int or None
            if idx is not None the index will be used for saving in the storage.
            This might overwrite already existing trajectories!

        Notes
        -----
        This also saves all contained frames in the snapshot if not done yet.
        A single Snapshot object can only be saved once!
        """
        self._set(idx, snapshot)

        reversed = snapshot._reversed
        snapshot._reversed = LoaderProxy(self, AbstractSnapshotStore.paired_idx(idx))
        reversed._reversed = LoaderProxy(self, idx)

        # mark reversed as stored
        self.index[reversed] = idx ^ 1

    def _init(self):
        """
        Initializes the associated storage to index configuration_indices in it
        """
        super(AbstractSnapshotStore, self)._init()

        self.init_variable('momentum_reversed', 'bool', chunksizes=(1,))

    def all(self):
        return Trajectory([LoaderProxy(self, idx) for idx in range(len(self))])


# =============================================================================================
# CONCRETE CLASSES FOR SNAPSHOT TYPES
# =============================================================================================

class SnapshotStore(AbstractSnapshotStore):
    """
    An ObjectStore for Snapshots in netCDF files.
    """

    def __init__(self):
        super(SnapshotStore, self).__init__(Snapshot)

    def to_dict(self):
        return {}

    def _set(self, idx, snapshot):
        self.vars['configuration'][idx] = snapshot.configuration
        self.vars['momentum'][idx] = snapshot.momentum
        self.write('configuration', idx ^ 1, snapshot)
        self.write('momentum', idx ^ 1, snapshot)

        super(SnapshotStore, self)._set(idx, snapshot)

    def _get(self, idx, from_reversed=False):
        if from_reversed:
            obj = self.cache[idx ^ 1]

            return Snapshot(
                configuration=obj.configuration,
                momentum=obj.momentum,
                is_reversed=not obj.is_reversed,
                topology=obj.topology,
                reversed_copy=LoaderProxy(self, AbstractSnapshotStore.paired_idx(idx))
            )
        else:
            configuration = self.vars['configuration'][idx]
            momentum = self.vars['momentum'][idx]
            momentum_reversed = self.vars['momentum_reversed'][idx]
            topology = self.storage.topology

            return Snapshot(
                configuration=configuration,
                momentum=momentum,
                is_reversed=momentum_reversed,
                topology=topology,
                reversed_copy=LoaderProxy(self, AbstractSnapshotStore.paired_idx(idx))
            )


    def _init(self):
        """
        Initializes the associated storage to index configuration_indices in it
        """
        super(SnapshotStore, self)._init()

        self.storage.create_store('configurations', ConfigurationStore())
        self.storage.create_store('momenta', MomentumStore())

        self.init_variable('configuration', 'lazyobj.configurations',
                           description="the snapshot index (0..n_configuration-1) of snapshot '{idx}'.",
                           chunksizes=(1,)
                           )

        self.init_variable('momentum', 'lazyobj.momenta',
                           description="the snapshot index (0..n_momentum-1) 'frame' of snapshot '{idx}'.",
                           chunksizes=(1,)
                           )

    # =============================================================================================
    # COLLECTIVE VARIABLE UTILITY FUNCTIONS
    # =============================================================================================

    @property
    def op_configuration_idx(self):
        """
        Returns aa function that returns for an object of this storage the idx

        Returns
        -------
        function
            the function that returns the idx of the configuration
        """

        def idx(obj):
            return self.index[obj.configuration]

        return idx

    @property
    def op_momentum_idx(self):
        """
        Returns aa function that returns for an object of this storage the idx

        Returns
        -------
        function
            the function that returns the idx of the configuration

        """

        def idx(obj):
            return self.index[obj.momentum]

        return idx

class ToySnapshotStore(AbstractSnapshotStore):
    """
    An ObjectStore for Snapshots in netCDF files.
    """

    def __init__(self):
        super(ToySnapshotStore, self).__init__(ToySnapshot)

    def to_dict(self):
        return {}

    def _set(self, idx, snapshot):
        self.vars['coordinates'][idx] = snapshot.coordinates
        self.vars['velocities'][idx] = snapshot.velocities
        self.write('coordinates', idx ^ 1, snapshot)
        self.write('velocities', idx ^ 1, snapshot)

        super(ToySnapshotStore, self)._set(idx, snapshot)

    def _get(self, idx, from_reversed=False):
        if from_reversed:
            obj = self.cache[idx ^ 1]

            return ToySnapshot(
                coordinates=obj.coordinates,
                velocities=obj.velocities,
                is_reversed=not obj.is_reversed,
                topology=obj.topology,
                reversed_copy=LoaderProxy(self, AbstractSnapshotStore.paired_idx(idx))
            )
        else:
            coordinates = self.vars['coordinates'][idx]
            velocities = self.vars['velocities'][idx]
            momentum_reversed = self.vars['momentum_reversed'][idx]

            return ToySnapshot(
                coordinates=coordinates,
                velocities=velocities,
                is_reversed=momentum_reversed,
                topology=self.storage.topology,
                reversed_copy=LoaderProxy(self, AbstractSnapshotStore.paired_idx(idx))
            )

    def _init(self):
        """
        Initializes the associated storage to index configuration_indices in it
        """
        super(ToySnapshotStore, self)._init()

        n_atoms = self.storage.n_atoms
        n_spatial = self.storage.n_spatial

        self.init_variable('coordinates', 'numpy.float32',
                           dimensions=('atom', 'spatial'),
                           description="coordinate of atom '{ix[1]}' in dimension " +
                                       "'{ix[2]}' of configuration '{ix[0]}'.",
                           chunksizes=(1, n_atoms, n_spatial)
                           )

        self.init_variable('velocities', 'numpy.float32',
                           dimensions=('atom', 'spatial'),
                           description="the velocity of atom 'atom' in dimension " +
                                       "'coordinate' of momentum 'momentum'.",
                           chunksizes=(1, n_atoms, n_spatial)
                           )



    # =============================================================================================
    # COLLECTIVE VARIABLE UTILITY FUNCTIONS
    # =============================================================================================

    @property
    def op_configuration_idx(self):
        """
        Returns aa function that returns for an object of this storage the idx

        Returns
        -------
        function
            the function that returns the idx of the configuration
        """

        def idx(obj):
            return self.index[obj.configuration]

        return idx

    @property
    def op_momentum_idx(self):
        """
        Returns aa function that returns for an object of this storage the idx

        Returns
        -------
        function
            the function that returns the idx of the configuration

        """

        def idx(obj):
            return self.index[obj.momentum]

        return idx


# =============================================================================================
# FEATURE BASED SINGLE CLASS FOR ALL SNAPSHOT TYPES
# =============================================================================================

class FeatureSnapshotStore(AbstractSnapshotStore):
    """
    An ObjectStore for Snapshots in netCDF files.
    """

    def __init__(self, snapshot_class):
        super(FeatureSnapshotStore, self).__init__(snapshot_class)

        self._variables = list()

        for feature in self.features:
            self._variables += getattr(ft, feature)._variables

    @property
    def features(self):
        return self.snapshot_class.__features__

    def _set(self, idx, snapshot):
        for variable in self._variables:
            self.vars[variable][idx] = getattr(snapshot, variable)
            self.write(variable, AbstractSnapshotStore.paired_idx(idx), snapshot)

        super(FeatureSnapshotStore, self)._set(idx, snapshot)

    def _get(self, idx, from_reversed=False):
        snapshot = self.snapshot_class.__new__(self.snapshot_class)
        if from_reversed:
            obj = self.cache[idx ^ 1]

            AbstractSnapshot.__init__(
                snapshot,
                not obj.is_reversed,
                LoaderProxy(
                    self,
                    AbstractSnapshotStore.paired_idx(idx)
                ),
                self.storage.topology
            )

            for variables in self._variables:
                setattr(snapshot, variables, getattr(obj, variables))

        else:
            AbstractSnapshot.__init__(
                snapshot,
                self.vars['momentum_reversed'][idx],
                LoaderProxy(
                    self,
                    AbstractSnapshotStore.paired_idx(idx)
                ),
                self.storage.topology
            )

            for variables in self._variables:
                setattr(snapshot, variables, self.vars[variables][idx])

        return snapshot

    def _init(self):
        super(FeatureSnapshotStore, self)._init()

        for feature in self.features:
            getattr(ft, feature)._init(self)
