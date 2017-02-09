# from openpathsampling.netcdfplus import UniqueNamedObjectStore
# from openpathsampling import CollectiveVariable
#
#
# class AttributeStore(UniqueNamedObjectStore):
#     """
#     ObjectStore to store a dict with StorableObject : value
#     """
#     def __init__(self):
#         super(CVStore, self).__init__(
#             CollectiveVariable
#         )
#
#     def _save(self, cv, idx):
#         self.vars['json'][idx] = cv
#
#         if cv.diskcache_enabled:
#             self.add_diskcache(cv)
#
#     def _load(self, idx):
#         op = self.vars['json'][idx]
#
#         cache_store = self.storage.snapshots.get_cv_cache(idx)
#
#         if cache_store is not None:
#             op.set_cache_store(cache_store)
#             op.diskcache_enabled = True
#             op.diskcache_chunksize = cache_store.chunksize
#             op.allow_incomplete = cache_store.allow_incomplete
#
#         return op
#
#     def sync(self, cv):
#         """
#         This will update the stored cache of the collective variable. It is
#         different from saving in that the object is only created if it is
#         saved (and the object caching will prevent additional creation)
#
#         Parameters
#         ----------
#         cv : :class:`openpathsampling.CollectiveVariable` or `None`
#             the objectdict to store. if `None` is given (default) then
#             all collective variables are synced
#
#         """
#         self.storage.snapshots.sync_cv(cv)
#
#     def complete(self, cv):
#         self.storage.snapshots.complete_cv(cv)
#
#     def sync_all(self):
#         map(self.sync, self)
#
#     def complete_all(self):
#         map(self.complete, self)
#
#     def add_diskcache(
#             self,
#             cv,
#             template=None,
#             allow_incomplete=None,
#             chunksize=None):
#         """
#         Return the storage.vars[''] variable that contains the values
#
#         Parameters
#         ----------
#         cv : :class:`openpathsampling.CollectiveVariable`
#             the objectdict you request the attached variable store
#         template : :obj:`openpathsampling.engines.BaseSnapshot`
#             an optional snapshot to be used to compute a test value of the CV.
#             This will determine the type and shape of the
#         allow_incomplete : bool
#             if `True` the added store can hold a part of all values. Useful if
#             the values are large and/or complex and you do not need them for all
#             snapshots
#         chunksize : int
#             for partial storage you can set a chunksize and speedup
#         """
#
#         if template is None:
#             if cv.diskcache_template is not None:
#                 template = cv.diskcache_template
#             elif len(self.storage.snapshots) > 0:
#                 template = self.storage.snapshots[0]
#
#             else:
#                 raise RuntimeError(
#                     'Need either at least one stored snapshot or a '
#                     'template snapshot to determine type and shape of the CV.')
#
#         self.storage.snapshots.add_cv(
#             cv,
#             template,
#             allow_incomplete=allow_incomplete,
#             chunksize=chunksize
#         )
#
#     def cache_store(self, cv):
#         """
#         Return the storage.vars[''] variable that contains the values
#
#         Parameters
#         ----------
#         cv : :class:`openpathsampling.CollectiveVariable`
#             the objectdict you request the attached variable store
#
#         Returns
#         -------
#         :class:`openpathsampling.netcdfplus.ObjectStore` or `netcdf4.Variable`
#
#         """
#         return self.storage.snapshots.cv_list[cv][0]
#
#     def has_cache(self, cv):
#         """
#         Test weather a CV has a diskstore attached
#
#         Parameters
#         ----------
#         cv : :obj:`openpathsampling.CollectiveVariable`
#             the CV you want to check
#
#         Returns
#         -------
#         bool
#             `True` if the CV has a diskstore attached
#         """
#         return cv in self.storage.snapshots.cv_list
#
#     def set_cache_store(self, cv):
#         """
#         Sets the attached storage of a CV to this store
#
#         Parameters
#         ----------
#         cv : :obj:`openpathsampling.CollectiveVariable`
#             the CV you want to set this store as its cache
#         """
#         if self.has_cache(cv):
#             cv.set_cache_store(self.cache_store(cv))
#         else:
#             raise RuntimeWarning(
#                 ('Your object is not stored in "%s" yet and hence a store ' +
#                  'for the cache cannot be attached.' +
#                  'Save your CV first and retry.') % self.storage)
#
#     def cache_all(self):
#         """
#         Fill the caches of all CVs
#
#         """
#
#         # load all CVs regularly
#         for cv in self:
#             # And cache the feature stores
#             store = self.storage.snapshots.cv_list.get(cv)
#             if store is not None:
#                 store[0].fill_cache()
