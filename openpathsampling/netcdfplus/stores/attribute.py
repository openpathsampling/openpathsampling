from .named import UniqueNamedObjectStore
from openpathsampling.netcdfplus.attribute import PseudoAttribute
from openpathsampling.netcdfplus.stores.value import ValueStore


class PseudoAttributeStore(UniqueNamedObjectStore):
    """
    ObjectStore to store additional attributes to an existing store
    """
    def __init__(self, obj_store=None):
        super(PseudoAttributeStore, self).__init__(
            PseudoAttribute
        )
        self.value_store = obj_store

    def _save(self, cv, idx):
        self.vars['json'][idx] = cv

        if cv.diskcache_enabled:
            self.add_diskcache(cv)
        else:
            self.vars['cache'][idx] = None

    def _load(self, idx):
        cv = self.vars['json'][idx]
        cache_store = None

        if 'cache' in self.vars:
            cache_store = self.vars['cache'][idx]

        if cache_store is not None:
            cv.set_cache_store(cache_store)
            cv.diskcache_enabled = True
            cv.diskcache_chunksize = cache_store.chunksize
            cv.allow_incomplete = cache_store.allow_incomplete

        return cv

    def key_store(self, cv):
        """
        Return the store that saves the key objects

        Parameters
        ----------
        cv : `PseudoAttribute`

        Returns
        -------
        `ObjectStore`
            the store the contains the key objects

        """
        return self.storage._objects[cv.key_class]

    def sync(self, cv):
        """
        This will update the stored cache of the attribute. It is
        different from saving in that the object is only created if it is
        saved (and the object caching will prevent additional creation)

        Parameters
        ----------
        cv : :class:`openpathsampling.netcdf.PseudoAttribute` or `None`
            the objectdict to store. if `None` is given (default) then
            all collective variables are synced

        """
        self.key_store(cv).sync_attribute(cv)

    def initialize(self):
        super(PseudoAttributeStore, self).initialize()
        self.create_variable('cache', 'obj.stores')

    def complete(self, cv):
        """
        Call the related complete function to precompute values for all objects

        Parameters
        ----------
        cv : `PseudoAttribute`

        """
        self.key_store(cv).complete_attribute(cv)

    def sync_all(self):
        """
        Sync all attribute caches to disk

        """
        map(self.sync, self)

    def complete_all(self):
        """
        Complete all attributes

        """
        map(self.complete, self)

    def add_diskcache(
            self,
            cv,
            template=None,
            allow_incomplete=None,
            chunksize=None):
        """
        Return the storage.vars[''] variable that contains the values

        Parameters
        ----------
        cv : :class:`openpathsampling.netcdf.PseudoAttribute`
            the objectdict you request the attached variable store
        template : :obj:`openpathsampling.engines.BaseSnapshot`
            an optional snapshot to be used to compute a test value of the CV.
            This will determine the type and shape of the
        allow_incomplete : bool
            if `True` the added store can hold a part of all values. Useful if
            the values are large and/or complex and you do not need them for all
            snapshots
        chunksize : int
            for partial storage you can set a chunksize and speedup
        """

        if template is None:
            if cv.diskcache_template is not None:
                template = cv.diskcache_template
            elif len(self.key_store(cv)) > 0:
                template = self.key_store(cv)[0]

            else:
                raise RuntimeError(
                    'Need either at least one stored object to use as a '
                    'template to determine type and shape of the CV. '
                    'No items in store: ' + str(self.key_store(cv)))

        self.key_store(cv).add_attribute(
            ValueStore,
            cv,
            template,
            allow_incomplete=allow_incomplete,
            chunksize=chunksize
        )

    def cache_store(self, cv):
        """
        Return the storage.vars[''] variable that contains the values

        Parameters
        ----------
        cv : :class:`openpathsampling.netcdf.PseudoAttribute`
            the objectdict you request the attached variable store

        Returns
        -------
        :class:`openpathsampling.netcdfplus.ObjectStore` or `netcdf4.Variable`

        """
        return self.key_store(cv).attribute_list[cv]

    def has_cache(self, cv):
        """
        Test weather a CV has a diskstore attached

        Parameters
        ----------
        cv : :obj:`openpathsampling.netcdf.PseudoAttribute`
            the CV you want to check

        Returns
        -------
        bool
            `True` if the CV has a diskstore attached
        """
        return cv in self.key_store(cv).attribute_list

    def set_cache_store(self, cv):
        """
        Sets the attached storage of a CV to this store

        Parameters
        ----------
        cv : :obj:`openpathsampling.netcdf.PseudoAttribute`
            the CV you want to set this store as its cache
        """
        if self.has_cache(cv):
            cv.set_cache_store(self.cache_store(cv))
        else:
            raise RuntimeWarning(
                ('Your object is not stored in "%s" yet and hence a store ' +
                 'for the cache cannot be attached.' +
                 'Save your CV first and retry.') % self.storage)

    def cache_all(self):
        """
        Fill the caches of all CVs

        """

        # load all CVs regularly
        for cv in self:
            # And cache the feature stores
            store = self.key_store(cv).attribute_list.get(cv)
            if store is not None:
                store.fill_cache()
