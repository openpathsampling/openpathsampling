import logging

import openpathsampling.engines as peng
from openpathsampling.netcdfplus import ObjectStore, \
    LRUChunkLoadingCache

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


class SnapshotValueStore(ObjectStore):
    def __init__(
            self,
            time_reversible=True,
            allow_incomplete=False,
            chunksize=256
    ):
        super(SnapshotValueStore, self).__init__(None)
        self.snapshot_index = None
        if not time_reversible and not allow_incomplete:
            raise RuntimeError(
                'Only time_reversible CVs can currently be '
                'stored using mode "complete"')

        self.time_reversible = time_reversible
        self.allow_incomplete = allow_incomplete
        self.chunksize = chunksize

        self.snapshot_pos = None
        self._len = 0

    def to_dict(self):
        return {
            'time_reversible': self.time_reversible,
            'allow_incomplete': self.allow_incomplete,
            'chunksize': self.chunksize
        }

    def create_uuid_index(self):
        return dict()

    def register(self, storage, prefix):
        super(SnapshotValueStore, self).register(storage, prefix)
        # print self.storage.__dict__.keys()
        self.snapshot_pos = self.storage.stores['snapshots'].pos

    def __len__(self):
        return len(self.variables['value'])

    # ==========================================================================
    # LOAD/SAVE DECORATORS FOR CACHE HANDLING
    # ==========================================================================

    def load(self, idx):
        pos = self.snapshot_pos(idx)
        # print idx.__uuid__ in self.storage.stores['snapshots'].index
        # print self.storage.stores['snapshots'].index[idx.__uuid__]
        # print self.snapshot_pos
        # print self.storage.stores['snapshots'].pos
        # print 'CV:', idx, pos

        if pos is None:
            return None

        if self.time_reversible:
            pos //= 2

        if self.allow_incomplete:
            # we want to load by uuid and it was not in cache.
            if pos in self.index:
                n_idx = self.index[pos]
            else:
                return None
            if n_idx < 0:
                return None
        else:
            if pos < self._len:
                n_idx = pos
            else:
                return None

        # if it is in the cache, return it
        try:
            obj = self.cache[n_idx]
            return obj

        except KeyError:
            pass

        obj = self.vars['value'][n_idx]

        self.cache[n_idx] = obj

        return obj

    def __setitem__(self, idx, value):
        pos = self.snapshot_pos(idx)

        if pos is None:
            return

        if self.time_reversible:
            pos //= 2

        if self.allow_incomplete:
            if pos in self.index:
                return

            n_idx = len(self.index)

            self.cache.update_size(n_idx)

        else:
            if pos < self._len:
                return

            n_idx = idx

        if self.allow_incomplete:
            # only if partial storage is used store index and update
            self.vars['index'][n_idx] = pos
            self.index[pos] = n_idx

        self.vars['value'][n_idx] = value
        self.cache[n_idx] = value
        self._len = max(self._len, n_idx + 1)

    def fill_cache(self):
        self.cache.load_max()

    def restore(self):
        if self.allow_incomplete:  # only if partial storage is used
            for pos, idx in enumerate(self.vars['index'][:]):
                self.index[idx] = pos

        self._len = len(self)
        self.initialize_cache()

    def initialize(self):
        self.initialize_cache()

    def initialize_cache(self):
        self.cache = LRUChunkLoadingCache(
            chunksize=self.chunksize,
            variable=self.vars['value']
        )
        self.cache.update_size()

    def __getitem__(self, item):
        # enable numpy style selection of objects in the store
        try:
            if isinstance(item, peng.BaseSnapshot):
                return self.load(item)
            elif type(item) is list:
                return [self.load(idx) for idx in item]
            elif item is Ellipsis:
                return iter(self)
        except KeyError:
            return None

    def get(self, item):
        if self.allow_incomplete:
            try:
                return self[item]
            except KeyError:
                return None
        else:
            return self[item]
