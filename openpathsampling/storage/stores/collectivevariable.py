from openpathsampling.netcdfplus import PseudoAttributeStore
import openpathsampling as paths


class CVStore(PseudoAttributeStore):
    """
    ObjectStore to store a dict with StorableObject : value
    """
    def __init__(self):
        super(CVStore, self).__init__()
        self.content_class = paths.BaseSnapshot

    def _load(self, idx):

        # repeated here to add backward compatibility for older CV support
        cv = self.vars['json'][idx]

        cache_store = None

        sn = self.storage.snapshots._get_cv_name(idx)
        if sn in self.storage.stores.name_idx:
            cache_store = self.storage.stores[sn]

        elif 'cache' in self.vars:
            cache_store = self.vars['cache'][idx]

        if cache_store is not None:
            cv.set_cache_store(cache_store)
            cv.diskcache_enabled = True
            cv.diskcache_chunksize = cache_store.chunksize
            cv.allow_incomplete = cache_store.allow_incomplete

        return cv
