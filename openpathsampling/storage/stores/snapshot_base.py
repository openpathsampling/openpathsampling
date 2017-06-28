import abc
import logging
from uuid import UUID

import openpathsampling.engines as peng
from openpathsampling.netcdfplus import IndexedObjectStore

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


# ==============================================================================
# ABSTRACT BASE CLASS FOR SNAPSHOTS
# ==============================================================================

class BaseSnapshotStore(IndexedObjectStore):
    """
    An ObjectStore for Snapshots in netCDF files.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, descriptor):
        """

        Attributes
        ----------
        descriptor : openpathsampling.engines.SnapshotDescriptor
            a descriptor knowing the snapshot class and a dictionary of
            dimensions and their lengths.

        """

        # Using a store with None as type will not interfere with the main
        # `SnapshotWrapperStore`
        super(BaseSnapshotStore, self).__init__(None, json=False)
        self.descriptor = descriptor
        self._dimensions = descriptor.dimensions
        self._cls = self.descriptor.snapshot_class

    @property
    def snapshot_class(self):
        return self._cls

    @property
    def dimensions(self):
        return self.descriptor['dimensions']

    def __repr__(self):
        return "store.%s[%s(%s)]" % (
            self.prefix, self._cls.__name__, 'BaseSnapshot')

    @staticmethod
    def paired_idx(idx):
        """
        Return the paired index

        Snapshots are stored in pairs (2n, 2n+1) where one is the reversed copy.
        This make storing CVs easier. This function allows to get the paired
        index or the index of snapshot.reversed

        The implementation uses the trick that all you have to do is flip the
        lowest bit that determines even or odd.

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

    def to_dict(self):
        return {
            'descriptor': self.descriptor,
        }

    def load(self, idx):
        pos = idx // 2

        # we want to load by uuid and it was not in cache.
        if pos in self.index:
            n_idx = self.index[pos]
        else:
            raise KeyError(idx)

        if n_idx < 0:
            return None

        obj = self._load(n_idx)

        if idx & 1:
            obj = obj.reversed

        return obj

    def _load(self, idx):
        """
        Load a snapshot from the storage.

        Parameters
        ----------
        idx : int
            the integer index of the snapshot to be loaded

        Returns
        -------
        snapshot : :obj:`openpathsampling.engines.BaseSnapshot`
            the loaded snapshot instance
        """

        # if not load and return it
        st_idx = int(idx)

        obj = self._cls.__new__(self._cls)
        self._cls.init_empty(obj)
        self._get(st_idx, obj)
        return obj

    def save(self, obj, idx=None):
        pos = idx // 2

        if pos in self.index:
            # has been saved so quit and do nothing
            return idx

        n_idx = len(self.index)

        # mark as saved so circular dependencies will not cause infinite loops
        self.index.append(pos)

        # make sure in nested saving that an IDX is not used twice!

        logger.debug('Saving ' + str(type(obj)) + ' using IDX #' + str(n_idx))

        try:
            self._save(obj, n_idx)
            self.vars['index'][n_idx] = pos

            # store the name in the cache
            self.cache[n_idx] = obj

        except:
            logger.debug('Problem saving %d !' % n_idx)
            # in case we did not succeed remove the mark as being saved
            del self.index[pos]
            raise

        self._set_id(n_idx, obj)

        return idx

    def _save(self, snapshot, idx):
        """
        Add the current state of the snapshot in the database.

        Parameters
        ----------
        snapshot :class:`openpathsampling.snapshots.AbstractSnapshot`
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

    @abc.abstractmethod
    def _get(self, idx, snapshot):
        pass

    @abc.abstractmethod
    def _set(self, idx, snapshot):
        pass

    def load_indices(self):
        self.index.extend(self.vars['index'])

    def all(self):
        return peng.Trajectory(map(self.proxy, self.index.list))

    def __iter__(self):
        for idx in range(0, len(self), 2):
            yield self[idx]

    def __getitem__(self, item):
        """
        Enable numpy style selection of object in the store
        """
        try:
            if type(item) is int or type(item) is str or type(item) is UUID:
                return self.load(item)
            elif type(item) is slice:
                return [self.load(idx)
                        for idx in range(*item.indices(len(self)))]
            elif type(item) is list:
                return [self.load(idx) for idx in item]
            elif item is Ellipsis:
                return iter(self)
        except KeyError:
            return None

    def __len__(self):
        return len(self.storage.dimensions[self.prefix]) * 2
