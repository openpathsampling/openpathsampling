from openpathsampling.netcdfplus import UniqueNamedObjectStore
from openpathsampling import CollectiveVariable


class CVStore(UniqueNamedObjectStore):
    """
    ObjectStore to store a dict with StorableObject : value
    """
    def __init__(self):
        super(CVStore, self).__init__(
            CollectiveVariable
        )

    def _save(self, cv, idx):
        self.vars['json'][idx] = cv

        if cv.diskcache_enabled:
            self.add_storage_caching(cv)

    def _load(self, idx):
        op = self.vars['json'][idx]
        cache_store = self.storage.snapshots.get_cv_cache(idx)
        if cache_store:
            op.set_cache_store(cache_store)

        return op

    def add_storage_caching(
            self,
            cv,
            template=None,
            allow_partial=None,
            chunksize=None):

        if template is None:
            if cv.diskcache_template is not None:
                template = cv.diskcache_template
            elif len(self.storage.snapshots) > 0:
                template = self.storage.snapshots[0]
            else:
                raise RuntimeError('Need either at least one stored snapshot or a '
                                   'template snapshot to determine type and shape of the CV.')

        self.storage.snapshots.add_cv(
            cv,
            template,
            allow_partial=allow_partial,
            chunksize=chunksize
        )

    def cache_store(self, cv):
        """
        Return the storage.vars[''] variable that contains the values

        Parameters
        ----------
        cv : :class:`openpathsampling.CollectiveVariable`
            the objectdict you request the attached variable store

        Returns
        -------
        :class:`openpathsampling.netcdfplus.ObjectStore` or `netcdf4.Variable`

        """
        return self.storage.snapshots.cv_list[cv][0]

    def has_cache(self, cv):
        return cv in self.storage.snapshots.cv_list

    def set_cache_store(self, cv):
        """
        Sets the attached storage of a CV to this store

        Parameters
        ----------
        cv : :obj:`openpathsampling.CollectiveVariable`
            the CV you want to set this store as its cache
        """
        if self.has_cache(cv):
            cv.set_cache_store(self.cache_store(cv))
        else:
            raise RuntimeWarning(('Your object is not stored in "%s" yet and hence a store ' +
                                  'for the cache cannot be attached.' +
                                 'Save your CV first and retry.') % self.storage)

    def sync(self, objectdict=None):
        """
        This will update the stored cache of the collective variable. It is
        different from saving in that the object is only created if it is
        saved (and the object caching will prevent additional creation)

        Parameters
        ----------
        objectdict : :class:`openpathsampling.CollectiveVariable` or `None`
            the objectdict to store. if `None` is given (default) then
            all collective variables are synced

        See also
        --------
        CollectiveVariable.sync

        """
        if objectdict is None:
            for cv in self:
                self.sync(cv)

            return

        objectdict.sync()

    def cache_all(self):
        """
        Fill the cache of all cvs

        """

        self.storage.snapshots.cache_all_cvs()
