__author__ = 'Jan-Hendrik Prinz'

import collections
import numpy as np

from openpathsampling.netcdfplus import LRUCache, LoaderProxy


class ChainDict(object):
    """
    A dictionary-like structure with a logic to generate values and pass unknown values to other instances

    The default ChainDict requires a list of keys. If you want to allow also single values you need to
    add a ChainDict that interprets single and iterables like ExpandSingle.

    The default for unknown keys is None. This is necessary to pass on what is unknown.
    Everything that is not known in the current instance is passed on to .post

    Examples
    --------
    Create a CachingDict
    >>> cd = CacheChainDict(LRUCache(2))
    >>> cd[[1, 2]]
    [None, None]

    There will be no result since there is no logic to generate values.
    >>> def f(x):
    >>>     print 'eval', x
    >>>     return x**2
    >>> fnc_cd = Function(f, requires_lists = False)
    >>> fnc_cd[[1, 2]]
    eval, 1
    eval, 2
    [1, 4]

    Combine both dicts to cache the results from the function
    >>> combo_cd = cd + fnc_cd
    >>> combo_cd[[1,2]]
    eval, 1
    eval, 2
    [1,4]

    First time the function is called and the cache filled. Second time the cache is used.

    >>> combo_cd[[1,2]]
    [1,4]

    Attributes
    ----------
    _post : ChainDict
        the ChainDict to be called when this instance cannot evaluate given keys

    """

    def __init__(self):
        self._post = None

    def __getitem__(self, items):
        # first apply the own _get functions to compute
        results = self._get_list(items)

        if self._post is not None:
            nones = [obj[0] for obj in zip(items, results) if obj[1] is None]
            if len(nones) == 0:
                return results
            else:
                rep = self._post[[p for p in nones]]
                self._set_list(nones, rep)

                it = iter(rep)
                return [it.next() if p[1] is None else p[1] for p in zip(items, results)]

        return results

    def __setitem__(self, key, value):
        if isinstance(key, collections.Iterable):
            self._set_list(key, value)
        else:
            if value is not None:
                self._set(key, value)

        # pass __setitem__ to underlying dicts as default
        if self._post is not None:
            self._post[key] = value

    def _set(self, item, value):
        """
        Implementation on how to set a single value to this chaindict

        Default implementation is to not store anything.
        This is mostly used in caching and stores
        """
        pass

    def _set_list(self, items, values):
        """
        Implementation on how to set multiple values to this chaindict

        Default implementation is to not store anything.
        This is mostly used in caching and stores
        """
        [self._set(item, value) for item, value in zip(items, values) if values is not None]

    def _get(self, item):
        """
        Implementation on how to get the value of a single key

        Default implementation returns None
        """
        # return None
        return None

    def _get_list(self, items):
        """
        Implementation on how to get the values of a list of keys

        Default is to use _get on each single key

        Returns
        -------
        list of object
        """
        return [self._get(item) for item in items]

    def __call__(self, items):
        return self[items]

    def __gt__(self, other):
        """
        Combine two ChainDicts first > next into a new one.

        Parameters
        ----------
        other : :class:`openpathsampling.chaindict.ChainDict`
            the chaindict to be attached as a fallback

        Returns
        -------
        :class:`openpathsampling.chaindict.ChainDict`
            the current object with the attached fallback chaindict.

        Examples
        --------
        >>> new_dict = first_dict > fall_back
        """
        last = self
        while last._post is not None:
            last = last._post

        last._post = other
        return self

    def __lt__(self, other):
        """
        Combine two ChainDicts seconds < first into a new one.

        Parameters
        ----------
        other : :class:`openpathsampling.chaindict.ChainDict`
            the chaindict to be attached as a fallback

        Returns
        -------
        :class:`openpathsampling.chaindict.ChainDict`
            the current object with the attached fallback chaindict.

        Examples
        --------

        >>> new_dict = fall_back < first_dict
        """
        last = other
        while last._post is not None:
            last = last.post

        last._post = self
        return other

    @property
    def passing_chain(self):
        """
        Return a list of chaindicts in order they will be tried.

        Returns
        -------
        list of :class:`openpathsampling.chaindict.ChainDict`
            the list of chaindicts in order they are called

        """
        chain = [self]
        while chain[-1]._post is not None:
            chain.append(chain[-1]._post)

        return chain

    def str_chain(self):
        """
        Return a string representation of the chain of dicts called.

        Returns
        -------
        str
            the string representation

        """
        return ' > '.join(map(lambda x : x.__class__.__name__, self.passing_chain))


class Wrap(ChainDict):
    """A ChainDict that passes on all requests to the underlying ChainDict

    """
    def __init__(self, post):
        """
        Parameters
        ----------
        post : :class:`openpathsampling.chaindict.ChainDict`
            the underlying chain dict to be used
        """
        super(Wrap, self).__init__()
        self._post = post

    def __getitem__(self, items):
        return self._post[items]

    def __setitem__(self, key, value):
        self._post[key] = value


class MergeNumpy(ChainDict):
    """All returned values from underlying ChainDicts will be turned into a numpy array
    """

    def __getitem__(self, items):
        return np.array(self._post[items])

    def __setitem__(self, key, value):
        self._post[key] = value


class ExpandSingle(ChainDict):
    """
    Iterables will be unrolled and passed as a list
    """

    def __getitem__(self, items):
        if type(items) is LoaderProxy:
            return self._post[[items]][0]
        if hasattr(items, '__iter__'):
            try:
                dummy = len(items)
            except AttributeError:
                # no length means unbound iterator and we cannot handle these
                raise AttributeError('Iterators that do not have __len__ implemented are not supported. ' +
                                'You can wrap your iterator in list() if you know that it will finish.')

            try:
                return self._post[items.as_proxies()]
            except AttributeError:
                # turn possible iterators into list since we have to do it anyway.
                # Iterators do not work
                return self._post[list(items)]

        else:
            return self._post[[items]][0]

    def __setitem__(self, key, value):
        self._post[key] = value


class Transform(ChainDict):
    """
    Applies a transformation to the input keys
    """
    def __init__(self, transform):
        """
        transform : function
            the function to be applied to all input keys on the underlying dicts
        """
        super(Transform, self).__init__()
        self.transform = transform

    def __getitem__(self, item):
        return self._post[self.transform(item)]

    def __setitem__(self, key, value):
        self._post[self.transform(key)] = value


class Function(ChainDict):
    """
    Uses a regular function to evaluate given keys.

    This works effective like a function called with square brackets
    """
    def __init__(self, fnc, requires_lists=True, scalarize_numpy_singletons=False):
        """
        Parameters
        ----------
        fnc : function
            the function to be evaluated to return values to keys
        requires_lists : bool
            if `True` we assume that it is faster to pass lists to this function instead
            of evaluating each key separately
        scalarize_numpy_singletons : bool
            if `True` eventual numpy objects that have length one in their last dimension, will
            be flattened by the last dimension. This is often useful if you have function that
            will by default return a list of results. In case your function does so, you can
            treat it as returning a scalar.

        """
        super(Function, self).__init__()
        self._eval = fnc
        self.requires_lists = requires_lists
        self.scalarize_numpy_singletons = scalarize_numpy_singletons

    def _get(self, item):
        if self._eval is None:
            return None

        if self.requires_lists:
            result = self._eval([item])[0]

        else:
            result = self._eval(item)

        if self.scalarize_numpy_singletons and result.shape[-1] == 1:
            return result.reshape(result.shape[:-1])

        return result

    def _get_list(self, items):
        if self._eval is None:
            return [None] * len(items)

        if self.requires_lists:
            results = self._eval(items)

            if self.scalarize_numpy_singletons and results.shape[-1] == 1:
                results = results.reshape(results.shape[:-1])

        else:
            results = [self._eval(obj) for obj in items]
            if self.scalarize_numpy_singletons and results[0].shape[-1] == 1:
                results = map(lambda x : x.reshape(x.shape[:-1]), results)

        return results

    def get_transformed_view(self, transform):
        def fnc(obj):
            return transform(self(obj))

        return fnc

    def __setitem__(self, key, value):
        # Cannot set the value of a function and it has no fallback
        pass


class CacheChainDict(ChainDict):
    """
    Return Values from a cache filled from returned values of the underlying ChainDicts and
    """
    def __init__(self, cache):
        """
        Parameters
        ----------
        cache : :class:`openpathsampling.netcdfplus.cache.Cache` or dict
            the cache to be used to store the data
        """
        super(CacheChainDict, self).__init__()
        self.cache = cache

    def _contains(self, item):
        return item in self.cache

    def _get(self, item):
        try:
            return self.cache[item]
        except KeyError:
            return None

    def _set(self, item, value):
        self.cache[item] = value


class ReversibleCacheChainDict(CacheChainDict):
    """
    Return Values from a cache filled from returned values of the underlying ChainDicts and
    """
    def __init__(self, cache, reversible=False):
        """
        Parameters
        ----------
        cache : :class:`openpathsampling.netcdfplus.cache.Cache` or dict
            the cache to be used to store the data
        """
        super(ReversibleCacheChainDict, self).__init__(cache)
        self.reversible = reversible

    def _set(self, item, value):
        self.cache[item] = value
        if self.reversible:
            self.cache[item.reversed] = value


class LRUChainDict(CacheChainDict):
    """
    Uses an LRUCache to cache values
    """
    def __init__(self, size_limit=1000000):
        """

        Parameters
        ----------
        size_limit : int
            the maximal allowed number of objects in the cache

        """
        super(LRUChainDict, self).__init__(LRUCache(size_limit))


class StoredDict(ChainDict):
    """
    ChainDict that has a store attached and returns existing values from the store
    """
    def __init__(self, key_store, value_store, main_cache, cache=None):
        """
        Parameters
        ----------
        key_store : storage.Store
            the store that references usable keys
        value_store : storage.Variable
            the store that references the store variable to store the values by index
        main_cache : :class:`openpathsampling.netcdfplus.cache.Cache` or dict
            the main cache used for non-stored objects so that the StoredDict can access
            values for these objects, too, if the objects has been saved in the meantime.
        cache : :class:`openpathsampling.netcdfplus.cache.Cache` or dict
            the cache used to access stored values faster. If `None` (default) an
            :class:`openpathsampling.netcdfplus.cache.LRUCache` with 1000000 (1M) entries is used.
        """
        super(StoredDict, self).__init__()
        self.value_store = value_store
        self.key_store = key_store
        self.main_cache = main_cache
        self.max_save_buffer_size = None
        if cache is None:
            cache = LRUCache(1000000)
        self.cache = cache
        self.storable = set()
        self._last_n_objects = 0

    def _set(self, item, value):
        key = self._get_key(item)
        if key is not None:
            self.cache[key] = value
            self.storable.add(key)

    def _set_list(self, items, values):
        [self._set(item, value) for item, value in zip(items, values) if value is not None]

        if self.max_save_buffer_size is not None and len(self.storable) > self.max_save_buffer_size:
            self.sync()

    def sync(self):
        # Sync objects that had been saved and afterwards the CV was computed
        if len(self.storable) > 0:
            keys = [idx for idx in sorted(list(self.storable)) if idx in self.cache]
            values = [self.cache[idx] for idx in keys]
            self.value_store[keys] = values
            self.storable.clear()

        # Sync objects that first had a value computed and were later stored
        # For these we need to check the main_cache

        if self._last_n_objects < len(self.key_store):
            keys = range(self._last_n_objects, len(self.key_store))
            objs = map(self.key_store.cache.get_silent, keys)
            values = map(self.main_cache.get_silent, objs)
            pairs = [(key, value) for key, value in zip(keys, values) if value is not None]
            keys, values = zip(*pairs)

            self.value_store[list(keys)] = list(values)
            for key, value in pairs:
                self.cache[key] = value
            self._last_n_objects = len(self.key_store)

    def cache_all(self):
        values = self.value_store[:]
        self.cache.clear()
        [self.cache.__setitem__(key, value) for key, value in enumerate(values)]

    def _get_key(self, item):
        if item is None:
            return None

        if hasattr(item, '_idx'):
            if item._store is self.key_store:
                return item._idx

        return self.key_store.index.get(item, None)

    def _get(self, item):
        key = self._get_key(item)

        if key is None:
            return None

        if key in self.cache:
            return self.cache[key]
        else:
            # update cache with specific strategy
            self.cache[key] = self.value_store[key]
            return self.cache[key]

    def _get_list(self, items):
        keys = map(self._get_key, items)

        idxs = [item for item in keys if item is not None and item not in self.cache]

        if len(idxs) > 0:
            sorted_idxs = sorted(list(set(idxs)))

            sorted_values = self.value_store[sorted_idxs]
            replace = [None if key is None else self.cache[key] if key in self.cache else
            sorted_values[sorted_idxs.index(key)] for key in keys]
        else:
            replace = [None if key is None else self.cache[key] if key in self.cache else None for key in keys]

        return replace


class ReversibleStoredDict(StoredDict):
    """
    ChainDict that has a store attached and return existing values from the store. Supports reversible items
    """
    def __init__(self, key_store, value_store, backward_store, main_cache, cache=None):
        """
        Parameters
        ----------
        key_store : storage.Store
            the store that references usable keys
        value_store : storage.Variable
            the store that references the store variable to store the values by index
        backward_store : storage.Variable
            the store that references the store variable to store the values by index
            only for reversed objects. If `backward_store` is `value_store` then the
            dict is assumed to be reversible in the sense that forward and backward
            objects have the same value
        main_cache : :class:`openpathsampling.netcdfplus.cache.Cache` or dict
            the main cache used for non-stored objects so that the StoredDict can access
            values for these objects, too, if the objects has been saved in the meantime.
        cache : :class:`openpathsampling.netcdfplus.cache.Cache` or dict
            the cache used to access stored values faster. If `None` (default) an
            :class:`openpathsampling.netcdfplus.cache.LRUCache` with 1000000 (1M) entries is used.

        """
        super(ReversibleStoredDict, self).__init__(key_store, value_store, main_cache)

        self._reversible = backward_store is value_store
        self.backward_store = backward_store

    def _set(self, item, value):
        key = self._get_key(item)
        if key is not None:
            self.cache[key] = value
            self.storable.add(key)
            if self._reversible:
                # if reversible store also for reversed
                s_key = key + 1 - 2 * (key % 2)
                self.cache[s_key] = value
                self.storable.add(s_key)

        if self.max_save_buffer_size is not None and len(self.storable) > self.max_save_buffer_size:
            self.sync()

    def sync(self):
        # Sync objects that had been saved and afterwards the CV was computed
        if len(self.storable) > 0:
            keys = [idx for idx in sorted(list(self.storable)) if idx in self.cache]
            keys_fw = [idx/2 for idx in keys if not idx & 1]
            keys_bw = [idx/2 for idx in keys if idx & 1]

            if keys_fw:
                values_fw = [self.cache[idx] for idx in keys if not idx & 1]
                self.value_store[keys_fw] = values_fw

            if keys_bw:
                values_bw = [self.cache[idx] for idx in keys if idx & 1]
                self.backward_store[keys_bw] = values_bw

            self.storable.clear()

        # Sync objects that first had a value computed and were later stored
        # For these we need to check the main_cache

        if self._last_n_objects < len(self.key_store):
            keys = range(self._last_n_objects, len(self.key_store))
            objs = map(self.key_store.cache.get_silent, keys)
            values = map(self.main_cache.get_silent, objs)

            if self._reversible:
                # double all pairs of values and remove Nones
                val_old = values
                values = map(lambda x : x[0] if x[0] is not None else x[1], zip(values[0::2], values[1::2]))
                values = [val for val in values for _ in (0, 1)]

            pairs = [(key, value) for key, value in zip(keys, values) if value is not None]
            if len(pairs) > 0:
                pair_fw = [(pair[0] / 2, pair[1]) for pair in pairs if not pair[0] & 1]
                if pair_fw:
                    keys_fw, values_fw = zip(*pair_fw)
                    self.value_store[list(keys_fw)] = list(values_fw)

                if not self._reversible:
                    pair_bw = [(pair[0] / 2, pair[1]) for pair in pairs if pair[0] & 1]
                    if pair_bw:
                        keys_bw, values_bw = zip(*pair_bw)
                        self.backward_store[list(keys_bw)] = list(values_bw)

                for key, value in pairs:
                    self.cache[key] = value

        self._last_n_objects = len(self.key_store)

    def cache_all(self):
        # TODO: This only makes sense if the cache can fit everything.

        values_fw = self.value_store[:]
        if self._reversible:
            values_bw = values_fw
        else:
            values_bw = self.backward_store[:]

        self.cache.clear()
        [self.cache.__setitem__(2*key, value) for key, value in enumerate(values_fw)]
        [self.cache.__setitem__(2*key + 1, value) for key, value in enumerate(values_bw)]

    def _get(self, item):
        key = self._get_key(item)

        if key is None:
            return None

        if key in self.cache:
            return self.cache[key]
        else:
            # update cache with specific strategy
            if not key & 1:
                val = self.value_store[key / 2]
            else:
                val = self.backward_store[key / 2]

            self.cache[key] = val

            if self._reversible:
                self.cache[key ^ 1] = val

            return val

    def _get_list(self, items):
        keys = map(self._get_key, items)

        idxs = [item for item in keys if item is not None and item not in self.cache]

        if len(idxs) > 0:
            sorted_idxs = sorted(list(set(idxs)))

            keys_fw = [int(idx/2) for idx in sorted_idxs if not idx & 1]
            keys_bw = [int(idx/2) for idx in sorted_idxs if idx & 1]

            if keys_fw:
                values_fw = self.value_store[keys_fw]
            if keys_bw:
                values_bw = self.backward_store[keys_bw]

            replace = [
                None if key is None else
                self.cache[key] if key in self.cache else
                values_fw[keys_fw.index(key/2)] if not key & 1 else
                values_bw[keys_bw.index(key/2)] for key in keys]
        else:
            replace = [None if key is None else self.cache[key] if key in self.cache else None for key in keys]

        return replace

