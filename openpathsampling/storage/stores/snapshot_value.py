import logging

import openpathsampling.engines as peng
from openpathsampling.netcdfplus import ValueStore

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


class SnapshotValueStore(ValueStore):
    def __init__(
            self,
            time_reversible=True,
            allow_incomplete=False,
            chunksize=256
    ):
        super(SnapshotValueStore, self).__init__(
            peng.BaseSnapshot,
            allow_incomplete=allow_incomplete,
            chunksize=chunksize)

        if not time_reversible and not allow_incomplete:
            raise RuntimeError(
                'Only time_reversible CVs can currently be '
                'stored using mode "complete"')

        self.time_reversible = time_reversible

    def to_dict(self):
        return {
            'time_reversible': self.time_reversible,
            'allow_incomplete': self.allow_incomplete,
            'chunksize': self.chunksize
        }

    def __len__(self):
        return len(self.variables['value'])

    # ==========================================================================
    # LOAD/SAVE DECORATORS FOR CACHE HANDLING
    # ==========================================================================

    def load(self, idx):
        pos = self.object_pos(idx)

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
        pos = self.object_pos(idx)

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
