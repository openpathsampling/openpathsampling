from __future__ import print_function
import openpathsampling as paths

from openpathsampling.netcdfplus import NamedObjectStore, ImmutableDictStore
import openpathsampling.engines as peng

from .stores import SnapshotWrapperStore
from .storage import Storage


class TrajectoryStorage(Storage):
    """
    A reduce storage that can only handle trajectories and CVs
    """

    def _create_storages(self):
        """
        Register all Stores used in the OpenPathSampling Storage

        """

        # objects with special storages

        self.create_store('trajectories', paths.storage.TrajectoryStore())
        self.create_store('snapshots', SnapshotWrapperStore())

        # normal objects

        self.create_store('topologies', NamedObjectStore(peng.Topology))
        self.create_store('engines', NamedObjectStore(peng.DynamicsEngine))

        # special stores

        self.create_store('tag', ImmutableDictStore())


class DistributedUUIDStorage(object):
    """
    A View on a storage that only changes the iteration over steps.

    Can be used for bootstrapping on subsets of steps and pass this object
    to analysis routines.

    """

    class MultiDelegate(object):
        """
        A delegate that will alter the ``iter()`` behaviour of the store

        Attributes
        ----------
        stores : list of :class:`openpathsampling.netcdfplus.ObjectStore`
            a list of references to the object stores used
        index : dict
            the dict containing all uuid references

        """

        def __init__(self, stores, index):
            self.stores = stores
            self.index = index

        def __iter__(self):
            for store, idx in self.index.values():
                yield store[idx]

        def __getitem__(self, item):
            store, idx = self.index.get(str(item), (None, None))
            if store is not None:
                return store[item]

        def __setitem__(self, key, value):
            pass

        def __len__(self):
            return len(self.index)

        def load(self, idx):
            return self[idx]

    def __init__(self, storages=None):
        """
        Parameters
        ----------

        storages : list of :class:`openpathsampling.storage.Storage`
            The storages to be used

        """
        self.storages = []

        self.main_storage = None
        self.stores = {}

        if storages is not None:
            map(self.add, storages)

    def add(self, storage):
        """
        Add a (read-only) file to the set of storages

        Parameters
        ----------
        storage

        """

        if type(storage) is list:
            [self.add(st) for st in storage]
            return

        if not hasattr(storage, 'use_uuid'):
            raise RuntimeError('The storage to be added does not use UUIDs!')

        for store in storage.stores:
            name = store.name
            if name not in self.stores:
                # a store is not yet present in the multi storage, so add it
                index = {}
                stores = []
                delegate = self.MultiDelegate(stores, index)
                setattr(self, name, delegate)
                self.stores[name] = delegate
            else:
                delegate = self.stores[str(name)]

            delegate.stores.append(store)

            for uuid, key in store.index.items():
                u = str(uuid)
                if u not in delegate.index:
                    delegate.index[u] = (store, key)

            # make all stores use the joint load functions
            store.register_fallback(getattr(self, name))
