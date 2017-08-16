import logging

from .object import ObjectStore
from openpathsampling.netcdfplus.cache import LRUChunkLoadingCache

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


class ValueStore(ObjectStore):
    """
    Store that stores a value by integer index

    Usually used to save additional attributes for objects

    See Also
    --------
    `PseudoAttribute`, `PseudoAttributeStore`
    """
    def __init__(
            self,
            key_class,
            allow_incomplete=False,
            chunksize=256
    ):
        super(ValueStore, self).__init__(None)
        self.key_class = key_class
        self.object_index = None
        self.allow_incomplete = allow_incomplete
        self.chunksize = chunksize

        self.object_pos = None
        self._len = 0

    def to_dict(self):
        return {
            'key_class': self.key_class,
            'allow_incomplete': self.allow_incomplete,
            'chunksize': self.chunksize
        }

    def create_uuid_index(self):
        return dict()

    def register(self, storage, prefix):
        super(ValueStore, self).register(storage, prefix)
        self.object_pos = self.storage._objects[self.key_class].pos

    def __len__(self):
        return len(self.variables['value'])

    # ==========================================================================
    # LOAD/SAVE DECORATORS FOR CACHE HANDLING
    # ==========================================================================

    def load(self, idx):
        pos = self.object_pos(idx)
        if pos is None:
            return None

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
        pos = self.object_pos(idx)

        if pos is None:
            return

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
            if isinstance(item, self.key_class):
                return self.load(item)
            elif type(item) is list:
                return [self.load(idx) for idx in item]
        except KeyError:
            pass

        return None

    def get(self, item):
        if self.allow_incomplete:
            try:
                return self.load(item)
            except KeyError:
                return None
        else:
            return self.load(item)
