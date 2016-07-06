from openpathsampling.netcdfplus import UniqueNamedObjectStore, ObjectStore
from openpathsampling.engines import BaseSnapshot
from openpathsampling import CollectiveVariable


class CVStore(UniqueNamedObjectStore):
    """
    ObjectStore to store a dict with StorableObject : value
    """
    def __init__(self):
        super(CVStore, self).__init__(
            CollectiveVariable
        )

    def initialize(self, units=None):
        super(CVStore, self).initialize()

        # index associated storage in class variable for all Trajectory instances to access

        if self.reference_by_uuid:
            self.storage.create_dimension('cv_cache')

            self.storage.create_variable(
                self.prefix + '_index',
                var_type='uuid',
                dimensions=('cv_cache',),
                description="the uuid for a specific cv_index",
                chunksizes=(10240,)
            )

    def to_dict(self):
        return {
        }

    def _save(self, objectdict, idx):
        """
        Save the current state of the cache to the storage.

        Parameters
        ----------
        objectdict : :class:`openpathsampling.CollectiveVariable`
            the objectdict to store
        idx : int
            the index
        """
        self.vars['json'][idx] = objectdict

    def add_storage_caching(
            self,
            cv,
            allow_partial=False,
            template=None):

        if template is None:
            if len(self.storage.snapshots) > 0:
                template = cv(self.storage.snapshots[0])
            else:
                raise RuntimeError('Need either at least one stored snapshot or a '
                                   'template snapshot to determine type and shapte of the CV.')

        self.storage.snapshots.add_cv(cv, cv(template), allow_partial=allow_partial)

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

    def _load(self, idx):
        op = self.vars['json'][idx]
        op.set_cache_store(self.storage.snapshots.get_cv_cache(idx))

        return op


class IntValueStore(ObjectStore):
    def __init__(self, cv):
        super(IntValueStore, self).__init__(None)

        self.cv = cv
        self._storable = None

    @property
    def uuid_ref(self):
        return self.storage.snapshots.index

    def create_int_index(self):
        return dict()

    # =============================================================================
    # LOAD/SAVE DECORATORS FOR CACHE HANDLING
    # =============================================================================

    def load(self, idx):
        """
        Returns an object from the storage.

        Parameters
        ----------
        idx : int
            the integer index of the object to be loaded

        Returns
        -------
        :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the loaded object
        """

        # we want to load by uuid and it was not in cache.
        if idx in self.index:
            n_idx = self.index[idx]
        else:
            if self.fallback_store is not None:
                return self.fallback_store.load(idx)
            elif self.storage.fallback is not None:
                return self.storage.fallback.stores[self.name].load(idx)
            else:
                raise ValueError('str %s not found in storage or fallback' % idx)

        if n_idx < 0:
            return None

        # if it is in the cache, return it
        try:
            obj = self.cache[n_idx]
            return obj

        except KeyError:
            pass

        obj = self.vars['value'][n_idx]

        self.index[idx] = n_idx
        self.cache[n_idx] = obj

        return obj

    def __setitem__(self, idx, value):
        """
        Saves an object to the storage.

        Parameters
        ----------
        idx : :py:class:`openpathsampling.engines.BaseSnapshot`
            the object to be stored
        value : anything that can be stored
            this includes storable objects, python numbers, numpy.arrays,
            strings, etc.

        """

        if idx in self.index:
            # has been saved so quit and do nothing
            return

        pos = self.storage.cvs.snapshot_index(idx)

        if pos is not None:

            n_idx = self.free()
            self.vars['value'][n_idx] = value
            self.vars['index'][n_idx] = pos

            self.index[idx] = n_idx
            self.cache[n_idx] = value

    def sync(self, cv):
        if not self.reference_by_uuid:
            # for uuids this cannot happen
            # necessary if we compute cvs that are not stored
            pass

    def restore(self):
        if self.reference_by_uuid:
            uuid_ref = self.uuid_ref
            for pos, idx in enumerate(self.vars['index'][:]):
                self.index[uuid_ref[idx]] = pos
        else:
            for pos, idx in enumerate(self.vars['index'][:]):
                self.index[idx] = pos

    def initialize(self):
        pass

    def __getitem__(self, item):
        """
        Enable numpy style selection of object in the store
        """
        try:
            if isinstance(item, BaseSnapshot):
                return self.load(item)
            elif type(item) is list:
                return [self.load(idx) for idx in item]
            elif item is Ellipsis:
                return self.iterator()
        except KeyError:
            return None

    def get(self, item):
        if item in self.index:
            return self[item]
        else:
            return None
