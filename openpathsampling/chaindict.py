__author__ = 'jan-hendrikprinz'

import collections


class ChainDict(dict):
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
        dict.__init__(self)
        self.post = None

    def __getitem__(self, items):
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
        self[items] = values

    def __setitem__(self, key, value):
        if isinstance(key, collections.Iterable):
            self._set_list(key, value)
        else:
            self._set(key, value)

    def _contains(self, item):
        return dict.__contains__(self, item)

    def _contains_list(self, items):
        return [dict.__contains__(self, item) for item in items]

    def _set(self, item, value):
        if value is not None:
            dict.__setitem__(self, item, value)

    def _set_list(self, items, values):
        [dict.__setitem__(self, item, value) for item, value in zip(items, values) if value is not None]

    def _get(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
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
        nones = values

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

class ExpandSingle(ChainDict):
    """
    Will take care of iterables
    """

    def __getitem__(self, items):
        if hasattr(items, '__iter__'):
            return self.post[items]
        else:
            return self.post[[items]][0]

    def __setitem__(self, key, value):
        self.post[key] = value

    def _add_new(self, items, values):
        pass

class ExpandMulti(ChainDict):
    """
    Will only request the unique keys to post
    """

    def __getitem__(self, items):
        if len(items) == 0:
            return []

        uniques = list(set(items))
        rep_unique = self.post[[p for p in uniques]]
        multi_cache = dict(zip(uniques, rep_unique))

        return [multi_cache[item] for item in items]

    def __setitem__(self, key, value):
        self.post[key] = value

    def _add_new(self, items, values):
        pass

class Transform(ChainDict):
    def __init__(self, transform):
        super(Transform, self).__init__()
        self.transform = transform

    def __getitem__(self, item):
        return self.post[self.transform(item)]

    def __setitem__(self, key, value):
        self.post[self.transform(key)] = value

    def _add_new(self, items, values):
        pass

class Function(ChainDict):
    def __init__(self, fnc, fnc_uses_lists=True):
        super(Function, self).__init__()
        self._fnc = fnc
        self.fnc_uses_lists = fnc_uses_lists

    def _eval(self, items):
        if hasattr(self, '_fnc'):
            return self._fnc(items)
        else:
            return None

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

class BufferedStore(Wrap):
    def __init__(self, name, dimensions, store, scope=None):
        self.storage = store.storage
        self._store = Store(name, dimensions, store, scope)
        self._cache = ChainDict()

        super(BufferedStore, self).__init__(
            post=self._store + self._cache
        )

    def sync(self):
        self._store.sync()

    def _add_new(self, items, values):
        for item, value in zip(items, values):
            if value is not None:
                if type(item) is tuple or len(item.idx) > 0 \
                        and self.storage in item.idx:
                    self._cache._set(item, value)
                    self._store._set(item, value)

class Store(ChainDict):
    def __init__(self, name, dimensions, store, scope=None):
        super(Store, self).__init__()
        self.name = name
        self.dimensions = dimensions
        self.store = store
        self.key_class = store.content_class

        if scope is None:
            self.scope = self
        else:
            self.scope = scope

        self.max_save_buffer_size = None

    def _add_new(self, items, values):
        [dict.__setitem__(self, item, value) for item, value in zip(items, values)]

        if self.max_save_buffer_size is not None and len(self) > self.max_save_buffer_size:
            self.sync()

    @property
    def storage(self):
        return self.store.storage

    def sync(self):
#        print 'Sync', len(self)
        storable = [ (key.idx[self.storage], value)
                            for key, value in self.iteritems()
                            if type(key) is not tuple and len(key.idx) > 0
                            and self.storage in key.idx]

        if len(storable) > 0:
            storable_sorted = sorted(storable, key=lambda x: x[0])
            storable_keys = [x[0] for x in storable_sorted]
            storable_values = [x[1] for x in storable_sorted]
            self.store.set_list_value(self.scope, storable_keys, storable_values)
            self.clear()
        else:
            self.clear()

    def _get_key(self, item):
        if item is None:
            return None

        if type(item) is tuple:
            if item[0].content_class is self.store.key_class:
                return item[1]
            else:
                return None

        elif self.storage in item.idx:
            return item.idx[self.storage]

        return None

    def _get(self, item):
        if dict.__contains__(self, item):
            return dict.__getitem__(self, item)

        key = self._get_key(item)

        if key is None:
            return None

        return self._load(key)

    def _get_list(self, items):
        cached, missing = self._split_list_dict(self, items)

        keys = [self._get_key(item) for item in missing]
        replace = self._apply_some_list(self._load_list, keys)

        return self._replace_none(cached, replace)

    def _load(self, key):
        return self.store.get_value(self.scope, key)

    def _load_list(self, keys):
        # This is to load all keys in ordered fashion since netCDF does not
        # allow reading in unsorted order using lists
        # TODO: Might consider moving this logic to the store, but this is faster
        # Also requesting an empty list raises an Error
        if len(keys) > 0:
            keys_sorted = sorted(enumerate(keys), key=lambda x: x[1])
            loadable_keys = [x[1] for x in keys_sorted]
            loadable_idxs = [x[0] for x in keys_sorted]
            values_sorted = self.store.get_list_value(self.scope, loadable_keys)
            ret = [0.0] * len(keys)
            [ret.__setitem__(idx, values_sorted[pos])
                    for pos, idx in enumerate(loadable_idxs)]
            return ret
        else:
            return []

    def _basetype(self, item):
        if type(item) is tuple:
            return item[0].content_class
        elif hasattr(item, 'base_cls'):
            return item.base_cls
        else:
            return type(item)

class MultiStore(Store):
    def __init__(self, store_name, name, dimensions, scope):
        super(Store, self).__init__()
        self.name = name
        self.dimensions = dimensions
        self.store_name = store_name
        self._storages = []

        if scope is None:
            self.scope = self
        else:
            self.scope = scope

        self.cod_stores = {}
        self.update_nod_stores()

    @property
    def storages(self):
        if hasattr(self.scope, 'idx'):
            if len(self.scope.idx) != len(self._storages):
                self._storages = self.scope.idx.keys()
            return self._storages
        else:
            return []

    def sync(self):
        if len(self.storages) != len(self.cod_stores):
            self.update_nod_stores()

        if len(self.cod_stores) == 0:
            return None

        [store.sync() for store in self.cod_stores.values()]

    def add_nod_store(self, storage):
        self.cod_stores[storage] = BufferedStore(
            self.name, self.dimensions, getattr(storage, self.store_name),
            self.scope
        )

    def update_nod_stores(self):
        for storage in self.cod_stores:
            if storage not in self.storages:
                del self.cod_stores[storage]

        for storage in self.storages:
            if storage not in self.cod_stores:
                self.add_nod_store(storage)

    def _add_new(self, items, values):
        if len(self.storages) != len(self.cod_stores):
            self.update_nod_stores()
        for s in self.cod_stores:
            self.cod_stores[s]._add_new(items, values)

    def _get(self, item):
        if len(self.storages) != len(self.cod_stores):
            self.update_nod_stores()

        if len(self.cod_stores) == 0:
            return None

        results = dict()
        for s in self.cod_stores:
            results[s] = self.cod_stores[s][item]

        output = None

        for result in results.values():
            if result is None:
                return None
            elif output is None:
                output = result

        return output


    def _get_list(self, items):
        if len(self.storages) != len(self.cod_stores):
            self.update_nod_stores()

        output = [None] * len(items)

        if len(self.cod_stores) == 0:
            return output

        results_list = dict()
        for s in self.cod_stores:
            results_list[s] = self.cod_stores[s][items]

        first = True
        for s, results in results_list.iteritems():
            if first:
                output = results
                first = False
            else:
                output = [None if item is None or result is None else output
                     for item, result in zip(output, results) ]

        return output



class UnwrapTuple(ChainDict):
    def __init__(self):
        super(UnwrapTuple, self).__init__()

    def __getitem__(self, items):
        return self.post([value[0].load(value[1])
            if type(value) is tuple else value for value in items])

    def __setitem__(self, key, value):
        self.post[key] = value