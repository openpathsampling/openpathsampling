__author__ = 'Jan-Hendrik Prinz'

import collections

import numpy as np
import weakref

from openpathsampling.storage.cache import LRUCache
from openpathsampling.storage.objproxy import LoaderProxy


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
    post : ChainDict
        the ChainDict to be called when this instance cannot evaluate given keys

    """

    def __init__(self):
        self.post = None

    def __getitem__(self, items):
        # first apply the own _get functions to compute
        results = self._get_list(items)

        if self.post is not None:
            nones = [obj[0] for obj in zip(items, results) if obj[1] is None]
            if len(nones) == 0:
                return results
            else:
                rep = self.post[[p for p in nones]]
                self._add_new(nones, rep)

                it = iter(rep)
                return [it.next() if p[1] is None else p[1] for p in zip(items, results)]

        return results

    def _add_new(self, items, values):
        pass

    def __setitem__(self, key, value):
        if isinstance(key, collections.Iterable):
            self._set_list(key, value)
        else:
            self._set(key, value)

    def _set(self, item, value):
        """
        Implementation on how to set a single value to this chaindict

        Default implementation is to not store anything.
        This is mostly used in caching and stores
        """
        pass

    def _set_list(self, items, values):
        """
        Implementation on how to set a list of keys to this chaindict

        Default is to call _set on all single keys
        """
        [self._set(item, value) for item, value in zip(items, values) if value is not None]

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
        """
        return [self._get(item) for item in items]

    def __call__(self, items):
        return self[items]

    ## ToDo: We might replace + by > to make the passing direction clear
    def __add__(self, other):
        """
        Combine two ChainDicts first + seconds into a new one.

        >>> new_dict = first_dict + fall_back
        """
        other.post = self
        return other


class Wrap(ChainDict):
    """A ChainDict that passes on all request to the underlying ChainDict
    """
    def __init__(self, post):
        """
        Parameters
        ----------
        post : chaindict
            the underlying chain dict to be used
        """
        super(Wrap, self).__init__()
        self.post = post

    def __getitem__(self, items):
        return self.post[items]

    def __setitem__(self, key, value):
        self.post[key] = value


class MergeNumpy(ChainDict):
    """All returned values from underlying ChainDicts will be turned into a numpy array
    """

    def __getitem__(self, items):
        return np.array(self.post[items])

class ExpandSingle(ChainDict):
    """
    Iterables will be unrolled and passed as a list
    """

    def __getitem__(self, items):
        if type(items) is LoaderProxy:
            return self.post[[items]][0]
        if hasattr(items, '__iter__'):
            try:
                dummy = len(items)
            except AttributeError:
                # no length means unbound iterator and we cannot handle these
                raise AttributeError('Iterators that do not have __len__ implemented are not supported. ' +
                                'You can wrap your iterator in list() if you know that it will finish.')

            try:
                return self.post[items.as_proxies()]
            except AttributeError:
                # turn possible iterators into list since we have to do it anyway.
                # Iterators do not work
                return self.post[list(items)]

        else:
            return self.post[[items]][0]

    def __setitem__(self, key, value):
        self.post[key] = value

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
        return self.post[self.transform(item)]

    def __setitem__(self, key, value):
        self.post[self.transform(key)] = value


class Function(ChainDict):
    """
    Uses a regular function to evaluate given keys.

    This works effective like a function called with square brackets
    """
    def __init__(self, fnc, requires_lists=True, single_as_scalar=False):
        """
        Parameters
        ----------
        fnc : function
            the function to be evaluated to return values to keys
        requires_lists : bool
            if true we assume that it is faster to pass lists to this function instead
            of evaluating each key separately
        """
        super(Function, self).__init__()
        self._eval = fnc
        self.requires_lists = requires_lists
        self.single_as_scalar = single_as_scalar

    def _contains(self, item):
        return False

    def _get(self, item):
        if self._eval is None:
            return None

        if self.requires_lists:
            result = self._eval([item])[0]

        else:
            result = self._eval(item)

        if self.single_as_scalar and result.shape[-1] == 1:
            return result.reshape(result.shape[:-1])

        return result

    def _get_list(self, items):
        if self._eval is None:
            return [None] * len(items)

        if self.requires_lists:
            results = self._eval(items)

            if self.single_as_scalar and results.shape[-1] == 1:
                results = results.reshape(results.shape[:-1])

        else:
            results = [self._eval(obj) for obj in items]
            if self.single_as_scalar and results[0].shape[-1] == 1:
                results = map(lambda x : x.reshape(x.shape[:-1]), results)

        return results


    def get_transformed_view(self, transform):
        def fnc(obj):
            return transform(self(obj))

        return fnc


class CacheChainDict(ChainDict):
    """
    Return Values from a cache filled from returned values of the underlying ChainDicts and
    """
    def __init__(self, cache):
        """
        Parameters
        ----------
        cache : dict-like class
            the dict or cache to be used to store the data
        """
        super(CacheChainDict, self).__init__()
        self.cache = cache

    def _contains(self, item):
        return item in self.cache

    def _set(self, item, value):
        if value is not None:
            self.cache[item] = value

    def _get(self, item):
        try:
            return self.cache[item]
        except KeyError:
            return None

    def _add_new(self, items, values):
        for item, value in zip(items, values):
            self.cache[item] = value
            if self.reversible:
                self.cache[item.reversed] = value

class ReversibleCacheChainDict(CacheChainDict):
    """
    Return Values from a cache filled from returned values of the underlying ChainDicts and
    """
    def __init__(self, cache, reversible=False):
        """
        Parameters
        ----------
        cache : dict-like class
            the dict or cache to be used to store the data
        """
        super(ReversibleCacheChainDict, self).__init__(cache)
        self.reversible = reversible

    def _add_new(self, items, values):
        for item, value in zip(items, values):
            self.cache[item] = value
            if self.reversible:
                self.cache[item.reversed] = value

class LRUChainDict(CacheChainDict):
    """
    Use a LRUCache to cache values
    """
    def __init__(self, size_limit=1000000):
        super(LRUChainDict, self).__init__(LRUCache(size_limit))


class StoredDict(ChainDict):
    """
    ChainDict that has a store attached and return existing values from the store
    """
    def __init__(self, key_store, value_store, main_cache, cache=None):
        """
        Parameters
        ----------
        key_store : storage.Store
            the store that references usable keys
        value_store : storage.Variable
            the store that references the store variable to store the values by index
        cache : dict-like or cache, default: None
            the cache used to access stored values faster. If None an LRUCache with
            100000 entries is used
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

    def _add_new(self, items, values):
        for item, value in zip(items, values):
            key = self._get_key(item)
            if key is not None:
                self.cache[key] = value
                self.storable.add(key)

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
    def __init__(self, key_store, value_store, main_cache, cache=None, reversible=True):
        """
        Parameters
        ----------
        key_store : storage.Store
            the store that references usable keys
        value_store : storage.Variable
            the store that references the store variable to store the values by index
        cache : dict-like or cache, default: None
            the cache used to access stored values faster. If None an LRUCache with
            100000 entries is used
        """
        super(ReversibleStoredDict, self).__init__(key_store, value_store, main_cache)
        self.reversible = reversible

    def _add_new(self, items, values):
        for item, value in zip(items, values):
            key = self._get_key(item)
            if key is not None:
                self.cache[key] = value
                self.storable.add(key)
                if self.reversible:
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
            values = [self.cache[idx] for idx in keys]
            self.value_store[keys] = values
            self.storable.clear()

        # Sync objects that first had a value computed and were later stored
        # For these we need to check the main_cache

        if self._last_n_objects < len(self.key_store):
            keys = range(self._last_n_objects, len(self.key_store))
            objs = map(self.key_store.cache.get_silent, keys)
            values = map(self.main_cache.get_silent, objs)

            if self.reversible:
                # double all pairs of values and remove Nones
                values = map(lambda x : x[0] or x[1], zip(values[0::2], values[1::2]))
                values = [val for val in values for _ in (0, 1)]

            pairs = [(key, value) for key, value in zip(keys, values) if value is not None]
            if len(pairs) > 0:
                keys, values = zip(*pairs)

                self._last_n_objects = len(self.key_store)

                self.value_store[list(keys)] = list(values)
                for key, value in pairs:
                    self.cache[key] = value
