__author__ = 'jan-hendrikprinz'

import collections

import numpy as np

from openpathsampling.storage.cache import LRUCache
from openpathsampling.storage.objproxy import LoaderProxy


class ChainDict(object):
    """
    Cache attached to Configuration indices stored in Configuration storage

    Parameters
    ----------
    name : string
        A short and unique name to be used in storage

    Attributes
    ----------
    name

    fnc : index (int) -> value (float)
        the function used to generate the cache values for the specified index. In essence a list
    dimensions : int
        the dimension of the stored values. Default is `1`
    content_class : Class
        the class of objects that can be stored
    fnc_uses_lists : boolean
        if True then in the case that the dict is called with several object at a time. The dict
        creates a list of missing ones and passes all of these to the evaluating function at once.
        Otherwise the fall-back is to call each item seperately. If possible always the multiple-
        option should be used.

    Attributes
    ----------
    content_class
    fnc_uses_lists
    dimensions
    """

    use_unique = True

    def __init__(self):
        super(ChainDict, self).__init__()
        self.post = None

    def __getitem__(self, items):
        results = self._get_list(items)

        # print 'Res', self.__class__.__name__, results

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

    def _contains(self, item):
        return dict.__contains__(self, item)

    def _contains_list(self, items):
        return [self._contains(item) for item in items]

    def _set(self, item, value):
        pass

    def _set_list(self, items, values):
        [self._set(item, value) for item, value in zip(items, values) if value is not None]

    def _get(self, item):
        return None

    def _get_list(self, items):
        return [ self._get(item) for item in items ]

    def __call__(self, items):
        return self[items]

    def __add__(self, other):
        other.post = self
        return other

    def _split_list_dict(self, dct, items):
        nones = [dct[item] if item in dct else None for item in items]
        missing = [item for item in items if item not in dct]

        return nones, missing

    def _split_list(self, keys, values):
        missing = [ obj[0] for obj in zip(keys,values) if obj[1] is None ]
        nones = [ None if obj[1] is None else obj[0] for obj in zip(keys,values)  ]

        return nones, missing

    def _apply_some_list(self, func, items):
        some = [item for item in items if item is not None]
        replace = func(some)
        it = iter(replace)

        return [ it.next() if obj is not None else None for obj in items ]

    def _replace_none(self, nones, replace):
        it = iter(replace)
        return [ obj if obj is not None else it.next() for obj in nones ]

class Wrap(ChainDict):
    def __init__(self, post):
        super(Wrap, self).__init__()
        self.post = post

    def __getitem__(self, items):
        return self.post[items]

    def __setitem__(self, key, value):
        self.post[key] = value

class MergeNumpy(ChainDict):
    """
    Will take care of iterables
    """

    def __getitem__(self, items):
        return np.array(self.post[items])

class ExpandSingle(ChainDict):
    """
    Will take care of iterables
    """

    def __getitem__(self, items):
        if type(items) is LoaderProxy:
            return self.post[[items]][0]
        if hasattr(items, '__iter__'):
            return self.post[items]
        else:
            return self.post[[items]][0]

    def __setitem__(self, key, value):
        self.post[key] = value

class ExpandMulti(ChainDict):
    """
    Will only request the unique keys to post
    """

    def __getitem__(self, items):
        return self.post[items]
        if len(items) == 0:
            return []


        uniques = list(set(items))
        rep_unique = self.post[[p for p in uniques]]
        multi_cache = dict(zip(uniques, rep_unique))

        return [multi_cache[item] for item in items]

    def __setitem__(self, key, value):
        self.post[key] = value

class Transform(ChainDict):
    def __init__(self, transform):
        super(Transform, self).__init__()
        self.transform = transform

    def __getitem__(self, item):
        return self.post[self.transform(item)]

    def __setitem__(self, key, value):
        self.post[self.transform(key)] = value

class Function(ChainDict):
    def __init__(self, fnc, fnc_uses_lists=True):
        super(Function, self).__init__()
        self._eval = fnc
        self.fnc_uses_lists = fnc_uses_lists

    def _contains(self, item):
        return False

    def _get(self, item):
        if self._eval is None:
            return None
#             raise KeyError('No cached values for item - %s' % str(item))

        if self.fnc_uses_lists:
            result = self._eval([item])
            return result[0]
        else:
            result = self._eval(item)
            return result

    def _get_list(self, items):
        if self._eval is None:
            return [None] * len(items)

        if self.fnc_uses_lists:
            result = self._eval(items)
            return result
        else:
            return [self._eval(obj) for obj in items]

    def get_transformed_view(self, transform):
        def fnc(obj):
            return transform(self(obj))

        return fnc

class CacheChainDict(ChainDict):
    def __init__(self, cache):
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


class LRUChainDict(CacheChainDict):
    def __init__(self, size_limit=1000000):
        super(LRUChainDict, self).__init__(LRUCache(size_limit))


class StoredDict(ChainDict):
    def __init__(self, key_store, value_store, cache=None):
        super(StoredDict, self).__init__()
        self.value_store = value_store
        self.key_store = key_store
        self.max_save_buffer_size = None
        if cache is None:
            cache = LRUCache(100000)
        self.cache = cache
        self.storable = set()

    def _add_new(self, items, values):
        for item, value in zip(items, values):
            key = self._get_key(item)
            if key is not None:
                self.cache[key] = value
                self.storable.add(key)

        if self.max_save_buffer_size is not None and len(self.storable) > self.max_save_buffer_size:
            self.sync()

    def sync(self):
        if len(self.storable) > 0:
            keys = [idx for idx in sorted(list(self.storable)) if idx in self.cache]
            values = [self.cache[idx] for idx in keys]
            self.value_store[keys] = values
            self.storable.clear()

    def cache_all(self):
        values = self.value_store[:]
        self.cache.clear()
        [self.cache.__setitem__(key, value) for key, value in enumerate(values)]

    def _get_key(self, item):
        if item is None:
            return None

        if hasattr(item, '_idx'):
            store, idx = item._idx.iteritems().next()
            if store is self.key_store:
                return idx

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
#        cached, missing = self._split_list(items, keys)

        idxs = [item for item in keys if item is not None and item not in self.cache]

        if len(idxs) > 0:
            sorted_idxs = sorted(list(set(idxs)))

            sorted_values = self.value_store[sorted_idxs]
            replace = [None if key is None else self.cache[key] if key in self.cache else
                            sorted_values[sorted_idxs.index(key)] for key in keys]
        else:
            replace = [None if key is None else self.cache[key] if key in self.cache else None for key in keys]

        return replace

class MultiStore(object):
    def __init__(self, store_name, name, dimensions, scope, unit=None):
        super(object, self).__init__()
        self.name = name
        self.dimensions = dimensions
        self.store_name = store_name
        self.unit = unit
        self._stores = []

        if scope is None:
            self.scope = self
        else:
            self.scope = scope

        self.cod_stores = {}
        self.update_nod_stores()

    @property
    def stores(self):
        return []
        if hasattr(self.scope, 'idx'):
            if len(self.scope.idx) != len(self._stores):
                self._stores = self.scope.idx.keys()
            return self._stores
        else:
            return []

    def sync(self):
        if len(self.stores) != len(self.cod_stores):
            self.update_nod_stores()

        if len(self.cod_stores) == 0:
            return None

        [store.sync() for store in self.cod_stores.values()]

    def cache_all(self):
        if len(self.stores) != len(self.cod_stores):
            self.update_nod_stores()

        if len(self.cod_stores) == 0:
            return None

        [store.cache_all() for store in self.cod_stores.values()]

    def add_nod_store(self, store):
        self.cod_stores[store] = StoredDict(
            self.name, self.dimensions, store,
            self.scope, unit=self.unit
        )

    def update_nod_stores(self):
        for store in self.cod_stores:
            if store not in self.stores:
                del self.cod_stores[store]

        for store in self.stores:
            if store not in self.cod_stores:
                self.add_nod_store(store)

    def _add_new(self, items, values):
        if len(self.stores) != len(self.cod_stores):
            self.update_nod_stores()
        for s in self.cod_stores:
            self.cod_stores[s]._add_new(items, values)

    def _get_list(self, items):
        if len(self.stores) != len(self.cod_stores):
            self.update_nod_stores()

        if len(self.cod_stores) == 0:
            return [None] * len(items)

        results_list = dict()
        for s in self.cod_stores:
            results_list[s] = self.cod_stores[s][items]

        first = True
        output = None
        for s, results in results_list.iteritems():
            if first:
                output = results
                first = False
            else:
                output = [None if item is None or result is None else item
                     for item, result in zip(output, results) ]

        return output