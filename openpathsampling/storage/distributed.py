import openpathsampling as paths
from openpathsampling.netcdfplus import NamedObjectStore, ImmutableDictStore
import openpathsampling.engines as peng

from openpathsampling.storage.storage import Storage


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
        self.create_store('snapshots', paths.storage.FeatureSnapshotStore(self._template.__class__))

        self.create_store('cvs', paths.storage.ReversibleObjectDictStore(
            paths.CollectiveVariable,
            peng.BaseSnapshot
        ))

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
        A delegate that will alter the ``iter()`` behaviour of the underlying store

        Attributes
        ----------
        stores : list of :class:`openpathsampling.netcdfplus.ObjectStore`
            a list of references to the object stores used
        uuid_idx : dict
            the dict containing all uuid references

        """

        def __init__(self, stores, uuid_idx):
            self.stores = stores
            self.uuid_idx = uuid_idx

        def __iter__(self):
            for store, idx in self.uuid_idx.itervalues():
                yield store[idx]

        def __getitem__(self, item):
            store, idx = self.uuid_idx.get(str(item), (None, None))
            if store is not None:
                return store[item]

        def __setitem__(self, key, value):
            pass

        def __len__(self):
            return len(self.uuid_idx)

        def load(self, idx):
            print idx
            return self[idx]

    def __init__(self, storages = None):
        """
        Parameters
        ----------

        storage : :class:`openpathsampling.storage.Storage`
            The storage the view is watching
        step_range : iterable
            An iterable object that species the step indices to be iterated over
            when using the view

        """
        self.storages = []

        self.template = None
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

        if self.template is None:
            self.template = storage.template
        else:
            if storage.template.__uuid__ != self.template.__uuid__:
                raise RuntimeWarning(
                    'The storage to be added uses a different template snapshot! It is NOT '
                    'recommended to join these storages'
                )

        for store in storage.stores:
            name = store.name
            if name not in self.stores:
                # a store is not yet present in the multi storage, so add it
                uuid_idx = {}
                stores = []
                delegate = self.MultiDelegate(stores, uuid_idx)
                setattr(self, name, delegate)
                self.stores[name] = delegate
            else:
                delegate = self.stores[str(name)]

            delegate.stores.append(store)

            for uuid, key in store.uuid_idx.iteritems():
                u = str(uuid)
                if u not in delegate.uuid_idx:
                    delegate.uuid_idx[u] = (store, key)

            # make all stores use the joint load functions
            store.register_fallback(getattr(self, name))
