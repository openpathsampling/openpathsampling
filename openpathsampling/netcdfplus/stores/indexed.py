from .object import ObjectStore

import logging

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


class IndexedObjectStore(ObjectStore):
    """
    ObjectStore storing objects in arbitrary order

    This has a prefilled .index which knows at which position a certain
    index is stored. This way you can circumvent holes and keep the file smaller
    """

    # ==========================================================================
    # LOAD/SAVE DECORATORS FOR CACHE HANDLING
    # ==========================================================================

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
            raise KeyError(idx)

        if n_idx < 0:
            return None

        # if it is in the cache, return it
        try:
            obj = self.cache[n_idx]
            return obj

        except KeyError:
            pass

        obj = self._load(n_idx)
        self.cache[n_idx] = obj

        return obj

    # def create_uuid_index(self):
    #     return dict()

    def save(self, obj, idx=None):
        """
        Saves an object to the storage.

        Parameters
        ----------
        obj : :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the object to be stored
        idx : int or string or `None`
            the index to be used for storing. This is highly discouraged since
            it changes an immutable object (at least in the storage). It is
            better to store also the new object and just ignore the
            previously stored one.

        """

        if idx in self.index:
            # has been saved so quit and do nothing
            return idx

        # n_idx = self.free()
        n_idx = len(self.index)

        # mark as saved so circular dependencies will not cause infinite loops
        self.index.append(idx)

        # make sure in nested saving that an IDX is not used twice!
        # self.reserve_idx(n_idx)

        logger.debug('Saving ' + str(type(obj)) + ' using IDX #' + str(n_idx))

        try:
            self._save(obj, n_idx)
            self.vars['index'][n_idx] = idx

            # store the name in the cache
            # if hasattr(self, 'cache'):
            self.cache[n_idx] = obj

        except:
            logger.debug('Problem saving %d !' % n_idx)
            # in case we did not succeed remove the mark as being saved
            del self.index[idx]
            # self.release_idx(n_idx)
            raise

        # self.release_idx(n_idx)
        self._set_id(n_idx, obj)

        return idx

    def restore(self):
        self.index.clear()
        self.index.extend(self.vars['index'][:])

    def initialize(self):
        super(IndexedObjectStore, self).initialize()

        self.create_variable(
            'index',
            'index'
        )
