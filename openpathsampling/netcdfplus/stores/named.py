from openpathsampling.netcdfplus.base import StorableNamedObject

from .object import ObjectStore

import logging

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')

import sys
if sys.version_info > (3, ):
    long = int


class NamedObjectStore(ObjectStore):
    def __init__(self, content_class, json=True, nestable=False):
        super(NamedObjectStore, self).__init__(
            content_class=content_class,
            json=json,
            nestable=nestable
        )

        self._names_loaded = False
        self._name_idx = dict()

        if self.content_class is not None \
                and not issubclass(self.content_class, StorableNamedObject):
            raise ValueError((
                'Content class "%s" must be subclassed from '
                'StorableNamedObject.') %
                self.content_class.__name__
            )

    def initialize(self):
        """
        Initialize the associated storage to allow for object storage. Mainly
        creates an index dimension with the name of the object.
        """
        super(NamedObjectStore, self).initialize()

        self.create_variable(
            "name", 'str',
            description='The name of the object',
            chunksizes=tuple([65536])
        )

    @property
    def name_idx(self):
        """
        Returns a dictionary of all names pointing to stored indices

        Returns
        -------
        dict of str : set
            A dictionary that has all stored names as keys and the values are a
            set of indices where an object with this name is found.

        """

        # if not done already cache names once
        # if not self._names_loaded:
        #     self.update_name_cache()

        return self._name_idx

    def load_indices(self):
        super(NamedObjectStore, self).load_indices()
        self.update_name_cache()

    def update_name_cache(self):
        """
        Update the internal name cache with all stored names in the store.

        This allows to load by name for named objects
        """
        if not self._names_loaded:
            for idx, name in enumerate(
                    self.storage.variables[self.prefix + "_name"][:]):
                self._update_name_in_cache(name, idx)

            self._names_loaded = True

    def _update_name_in_cache(self, name, idx):
        # make sure to cast unicode to str
        name = str(name)
        if name:
            if name not in self._name_idx:
                self._name_idx[name] = {idx}
            else:
                if idx not in self._name_idx[name]:
                    self._name_idx[name].add(idx)

    def cache_all(self):
        """Load all samples as fast as possible into the cache

        """
        if not self._cached_all:
            idxs = range(len(self))
            jsons = self.variables['json'][:]
            names = self.variables['name'][:]

            map(self.add_single_to_cache_named,
                idxs,
                names,
                jsons)

            # [self.add_single_to_cache(i, n, j) for i, n, j in zip(
            #     idxs,
            #     names,
            #     jsons)]

            self._cached_all = True

    def add_single_to_cache_named(self, idx, name, json):
        """
        Add a single object to cache by json

        Parameters
        ----------
        idx : int
            the index where the object was stored
        name : str
            the name of the object if it exists
        json : str
            json string the represents a serialized version of the stored object
        """

        if idx not in self.cache:
            obj = self.simplifier.from_json(json)

            self._get_id(idx, obj)

            self.cache[idx] = obj
            self.index[obj.__uuid__] = idx

            setattr(obj, '_name', name)

            return obj

    def find(self, name):
        """
        Return last object with a given name

        Parameters
        ----------
        name : str
            the name to be searched for

        Returns
        -------
        :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the last object with a given name. This is to mimic immutable
            objects. Once you (re-)save with the same name you replace the
            old one and hence you leed to load the last stored one.

        """
        return self.load(name)

    def find_indices(self, name):
        """
        Return indices for all objects with a given name

        Parameters
        ----------
        name : str
            the name to be searched for

        Returns
        -------
        list of int
            a list of indices in the storage for all found objects,
            can be empty [] if no objects with that name exist

        """
        return sorted(list(self._name_idx[name]))

    def find_all(self, name):
        if len(self._name_idx[name]) > 0:
            return self[sorted(list(self._name_idx[name]))]

    # ==========================================================================
    # LOAD/SAVE DECORATORS FOR CACHE HANDLING
    # ==========================================================================

    def load(self, idx):
        """
        Returns an object from the storage.

        Parameters
        ----------
        idx : int or str
            either the integer index of the object to be loaded or a string
            (name) for named objects. This will always return the last object
            found with the specified name. This allows to effectively change
            existing objects.

        Returns
        -------
        :py:class:`openpathsampling.netcdfplus.base.StorableNamedObject`
            the loaded object
        """

        if type(idx) is str:
            # we want to load by name and it was not in cache.
            if idx in self.name_idx:
                if len(self.name_idx[idx]) > 1:
                    if self._log_debug:
                        logger.debug((
                            'Found name "%s" multiple (%d) times in storage! '
                            'Loading last!') % (
                            idx, len(self.cache[idx])))

                n_idx = sorted(list(self.name_idx[idx]))[-1]
            else:
                raise ValueError('str "' + idx + '" not found in storage')

        # --- start super of ObjectStore ---

        elif isinstance(idx, (int, long)):
            if idx < 1000000000:
                n_idx = idx
            elif idx in self.index:
                n_idx = self.index[idx]
            else:
                if self.fallback_store is not None:
                    return self.fallback_store.load(idx)
                elif self.storage.fallback is not None:
                    return self.storage.fallback.stores[self.name].load(idx)
                else:
                    raise ValueError(
                        'str %s not found in storage or fallback' % idx)

        # elif type(idx) is not int:
        #     raise ValueError((
        #          'indices of type "%s" are not allowed in named storage '
        #          '(only str and int)') % type(idx).__name__
        #                      )
        # else:
        #     n_idx = int(idx)

        if n_idx < 0:
            return None

        # if it is in the cache, return it
        try:
            obj = self.cache[n_idx]
            if self._log_debug:
                logger.debug(
                    'Found IDX #' + str(idx) + ' in cache. Not loading!')
            return obj

        except KeyError:
            pass

        if self._log_debug:
            logger.debug(
                'Calling load object of type `%s` @ IDX #%d' %
                (self.content_class.__name__, n_idx))

        if n_idx >= len(self):
            logger.warning(
                'Trying to load from IDX #%d > number of object %d' %
                (n_idx, len(self)))
            return None
        elif n_idx < 0:
            logger.warning((
                               'Trying to load negative IDX #%d < 0. '
                               'This should never happen!!!') % n_idx)
            raise RuntimeError(
                'Loading of negative int should result in no object. '
                'This should never happen!')
        else:
            obj = self._load(n_idx)

        if self._log_debug:
            logger.debug(
                'Calling load object of type %s and IDX # %d ... DONE' %
                (self.content_class.__name__, n_idx))

        if obj is not None:
            self._get_id(n_idx, obj)

            setattr(obj, '_name',
                    self.storage.variables[self.prefix + '_name'][n_idx])
            # make sure that you cannot change the name of loaded objects
            obj.fix_name()

            # finally store the name of a named object in cache
            # self._update_name_in_cache(obj._name, n_idx)

            # update cache there might have been a change due to naming
            self.cache[n_idx] = obj

            if self._log_debug:
                logger.debug(
                    'Try loading UUID object of type %s and IDX # %d ... DONE' %
                    (self.content_class.__name__, n_idx))

        if self._log_debug:
            logger.debug(
                'Finished load object of type %s and IDX # %d ... DONE' %
                (self.content_class.__name__, n_idx))

        # --- end ---

        return obj

    def save(self, obj, idx=None):
        """
        Saves an object to the storage.

        Parameters
        ----------
        obj : :py:class:`openpathsampling.netcdfplus.base.StorableNamedObject`
            the object to be stored
        idx : int or string or `None`
            the index to be used for storing. This is highly discouraged since
            it changes an immutable object (at least in the storage). It is
            better to store also the new object and just ignore the
            previously stored one.

        """

        is_str = type(idx) is str

        if not is_str and idx is not None:
            raise ValueError(
                'Unsupported index type (only str or None allowed).')

        name = obj._name

        if is_str:
            obj.name = idx
            name = obj._name

        if name is None:
            # this should not happen!
            logger.debug(
                "Nameable object has not been initialized correctly. "
                "Has None in _name")
            raise AttributeError(
                '_name needs to be a string for nameable objects.')

        obj_fixed = obj._name_fixed
        obj_name = obj._name
        # we fix the name just in case we try in recursive saving to store
        # one object twice with different names. Storing with the same name
        # is fine!
        obj.fix_name()

        try:
            reference = super(NamedObjectStore, self).save(obj)
        except:
            # if saving did not work unlock the name if is was un-fixed before
            obj._name_fixed = obj_fixed
            obj._name = obj_name
            raise

        n_idx = self.index[obj.__uuid__]
        self.storage.variables[self.prefix + '_name'][n_idx] = name
        self._update_name_in_cache(name, n_idx)

        return reference


class UniqueNamedObjectStore(NamedObjectStore):

    # ==========================================================================
    # LOAD/SAVE DECORATORS FOR CACHE HANDLING
    # ==========================================================================

    def __init__(self, content_class, json=True, nestable=False):
        super(UniqueNamedObjectStore, self).__init__(
            content_class=content_class,
            json=json,
            nestable=nestable)

        self._free_name = set()

    def reserve_name(self, name):
        """
        Locks a name as used

        Parameters
        ----------
        name : str
            the name to be locked for storage
        """
        if name != "":
            self._free_name.add(name)

    def release_name(self, name):
        """
        Releases a locked name

        Parameters
        ----------
        name : str
            the name to be released for being used as a name
        """
        self._free_name.discard(name)

    def is_name_locked(self, name):
        """
        Test whether in a unique name store a name is already taken

        Parameters
        ----------
        name : str or `None`
            the name to be tested.

        Returns
        -------
        bool
            the result of the test. If the name exists or is reserved during
            a saving event this will return `True` and return `False` if the
            name is free.

        """
        if name is None:
            return False

        return name in self.name_idx or name in self._free_name

    def save(self, obj, idx=None):
        """
        Saves an object to the storage.

        Parameters
        ----------
        obj : :py:class:`openpathsampling.netcdfplus.base.StorableNamedObject`
            the object to be stored
        idx : string or `None`
            the index to be used for storing. This is highly discouraged since
            it changes an immutable object (at least in the storage). It is
            better to store also the new object and just ignore the
            previously stored one.

        """

        is_str = type(idx) is str

        if not is_str and idx is not None:
            raise ValueError(
                'Unsupported index type (only str or None allowed).')

        name = obj._name
        fixed = obj._name_fixed
        err = list()

        if is_str:
            if fixed:
                if name != idx:
                    # saving fixed under different name is not possible.
                    if obj in self:
                        err.append(
                            ('Cannot rename object to "%s". '
                             'Already saved with name "%s" !') % (idx, name)
                        )
                    else:
                        err.append(
                            ('Cannot rename object to "%s". '
                             'Already fixed name "%s" !') % (idx, name)
                        )

                        if self.is_name_locked(name):
                            err.append(
                                ('Current name "%s" is also already taken in '
                                 'unique name store. This means you cannot '
                                 'save object "%s" at all. '
                                 'In general this should not happen to unsaved '
                                 'objects unless you fixed the name of the '
                                 'object yourself. Check your code '
                                 'for the generation of objects of the same '
                                 'name.') %
                                (name, obj)
                            )
                        else:
                            err.append(
                                ('Current name "%s" is still free. Saving '
                                 'without giving a specific name '
                                 'should work. If that is what you want '
                                 'to do.') % name
                            )
                else:
                    # already fixed, but with same name. Okay. Check if stored
                    if obj in self.index:
                        return self.reference(obj)
            else:
                # name is not fixed yet. Check, if we can save or name is taken
                if self.is_name_locked(idx):
                    err.append(
                        ('New name "%s" already taken in unique name store. ' +
                         'Try different name instead.') % idx
                    )

                    if self.is_name_locked(name):
                        err.append((
                            'Current name "%s" already taken in unique '
                            'name store. ') % name
                        )
                    else:
                        err.append(
                            ('Current name "%s" is still free. Saving without '
                             'giving a specific name should work') % name
                        )
        else:
            if fixed:
                # no new name, but fixed. Check if already stored.
                if obj.__uuid__ in self.index:
                    return self.reference(obj)

                # if not stored yet check if we could
                if self.is_name_locked(name):
                    err.append(
                        ('Current name "%s" is already taken in unique name '
                         'store. This means you cannot save object "%s" at '
                         'all. In general this should not happen to unsaved '
                         'objects unless you fixed the name of the object '
                         'yourself. Check your code for the generation of '
                         'objects of the same name.') %
                        (name, obj)
                    )
            else:
                # no new name and not fixed. Just check if current name is taken
                if self.is_name_locked(name):
                    err.append(
                        ('Current name "%s" is already taken in unique name '
                         'store %s. Try renaming object or saving using other '
                         'name.') % (name, self.name)
                    )

        if len(err) > 0:
            raise RuntimeWarning('/n'.join(err))

        # no errors, reserve the name for nested saving and actually call save
        self.reserve_name(name)

        try:
            reference = super(UniqueNamedObjectStore, self).save(obj, idx)
        finally:
            self.release_name(name)

        return reference
