from .named import NamedObjectStore

from future.utils import iterkeys

import logging

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


class DictStore(NamedObjectStore):
    def __init__(self):
        super(DictStore, self).__init__(
            None,
            json='json'
        )

    def to_dict(self):
        return {}

    def load(self, idx):
        """
        Returns an object from the storage.

        Parameters
        ----------
        idx : str
            a string (name) of the objects. This will always return the
            last object found with the specified name. If immutable is true
            for the store it assures that there is only a single object per name

        Returns
        -------
        :class:`openpathsampling.netcdfplus.base.StorableObject`
            the loaded object
        """

        if type(idx) is str:
            n_idx = -1

            # we want to load by name and it was not in cache.
            if idx not in self.name_idx:
                logger.debug('Name "%s" not found in the storage!' % idx)
                raise KeyError('str "' + idx + '" not found in storage')

            if idx in self.name_idx:
                if len(self.name_idx[idx]) > 1:
                    logger.debug((
                        'Found name "%s" multiple (%d) times in storage! '
                        'Loading last!') % (
                        idx, len(self.name_idx[idx])))

                n_idx = sorted(list(self.name_idx[idx]))[-1]

                logger.debug(
                    'Found name "%s" in storage! Loading IDX %d!' %
                    (idx, n_idx))

        elif type(idx) is int:
            n_idx = idx
        else:
            raise ValueError(
                'Unsupported index type (only str and int allowed).')

        # turn into python int if it was a numpy int (in some rare cases!)
        n_idx = int(n_idx)

        logger.debug(
            'Calling load object of type `%s` and @ IDX #%d' %
            (str(self.content_class), n_idx))

        if n_idx >= len(self):
            logger.warning(
                'Trying to load from IDX #' + str(n_idx) +
                ' > number of objects ' + str(len(self))
            )
            raise RuntimeError(
                'Loading of too large int should be attempted. '
                'Problem in name cache. This should never happen!')
        elif n_idx < 0:
            logger.warning(
                'Trying to load negative IDX #' + str(n_idx) + ' < 0. '
                'This should never happen!!!'
            )
            raise RuntimeError(
                'Loading of negative int should result in no object. '
                'This should never happen!'
            )
        else:

            logger.debug('Loading named object from index IDX # %d' % n_idx)
            obj = self._load(n_idx)

            logger.debug(
                'Loading named object from index IDX # %d.. DONE' % n_idx)

        return obj

    def restore(self):
        self.update_name_cache()

    def save(self, obj, idx=None):
        """
        Saves an object to the storage.

        Parameters
        ----------
        obj : :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the object to be stored
        idx : string or `None`
            the string index to be used for storing. Objects will not be
            replaced but stored again with the same name. When loading the
            last stored object under the idx is retrieved. Effectively mimicking
            a mutual dict with versioning. We usually encourage for most cases
            to use the immutual dict class
            :class:`openpathsampling.netcdfplus.ImmutableDictStore` instead
            to avoid ambiguity in stored objects.

        See Also
        --------
        :class:`openpathsampling.netcdfplus.ImmutableDictStore`

        """

        if idx is None:
            # a DictStore needs a specific name
            raise ValueError(
                'Saving in a DictStore without specifying a '
                'string key is not allowed. ')

        if type(idx) is not str:
            # key needs to be a string
            raise ValueError(
                'Index `%s` for DictStore needs to be a string! ' % idx)

        n_idx = int(self.free())
        # make sure in nested saving that an IDX is not used twice!
        # self.reserve_idx(n_idx)

        logger.debug(
            'Saving `%s` with name `%s` @ IDX #%d' %
            (str(obj.__class__), idx, n_idx))
        self._save(obj, n_idx)

        self.storage.variables[self.prefix + '_name'][n_idx] = idx
        self._update_name_in_cache(idx, n_idx)

        return n_idx

    def keys(self):
        return self.name_idx.keys()

    # def iterkeys(self):
    #     return self.name_idx.keys()

    def __iter__(self):
        return iter(iterkeys(self.name_idx))

    def iteritems(self):
        for name in self:
            yield name, self[name]

    def get(self, idx, default=None):
        try:
            return self.load(idx)
        except KeyError:
            return default


class ImmutableDictStore(DictStore):

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

        if idx in self.name_idx:
            # immutable means no duplicates, so quit
            raise RuntimeWarning(
                'Cannot re-save existing key "%s" in '
                'immutable dict store.' % idx
            )

        return super(ImmutableDictStore, self).save(obj, idx)
