__author__ = 'jan-hendrikprinz'

import collections


class NestableObjectDict(dict):
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

    def __iter__(self):
        return None

    def __getitem__(self, items):
        is_listable = isinstance(items, collections.Iterable)

        if is_listable:
            results = self._get_list(items)
        else:
            results = self._get(items)

        if self.post is not None:
            if is_listable:
                nones = [obj[0] for obj in zip(items, results) if obj[1] is None]
                if len(nones) == 0:
                    return results
                else:
                    rep = self.post[[p for p in nones]]
                    self._add_new(nones, rep)

                    it = iter(rep)
                    return [it.next() if p[1] is None else p[1] for p in zip(items, results)]
            else:
                if results is None:
                    rep = self.post[items]

                    self._add_new(items, rep)
                    return rep
                else:
                    return results

        return results

    def _add_new(self, items, values):
        self[items] = values

    def __setitem__(self, key, value):
        if isinstance(key, collections.Iterable):
            self._set_list(key, value)
        else:
            self._set(key, value)

#    def __contains__(self, item):
#        return dict.__contains__(self, item) or self.in_store(item)


    def _contains(self, item):
        return dict.__contains__(self, item)

    def _contains_list(self, items):
        return [dict.__contains__(self, item) for item in items]

    def _set(self, item, value):
        dict.__setitem__(self, item, value)

    def _set_list(self, items, values):
        [dict.__setitem__(self, item, value) for item, value in zip(items, values)]

    def _get(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            return None

    def _get_list(self, items):
        return [ self._get(item) for item in items ]

    def push(self):
        pass


    def existing(self, objs):
        """
        Find a subset of indices that are present in the cache

        Parameters
        ----------
        indices : list of int
            the initial list of indices to be tested

        Returns
        -------
        existing : list of int
            the subset of indices present in the cache
        """
        return [obj for obj in objs if obj in self]

    def missing(self, objs):
        """
        Find a subset of indices that are NOT present in the cache

        Parameters
        ----------
        indices : list of int
            the initial list of indices to be tested

        Returns
        -------
        existing : list of int
            the subset of indices NOT present in the cache
        """
        return [obj for obj in objs if obj not in self]

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


class CODExpandMulti(NestableObjectDict):
    """
    Will only request the unique keys to post
    """

    def __getitem__(self, items):
        is_list = isinstance(items, collections.Iterable)

        if not is_list:
            return self.post[items]

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

class CODTransform(NestableObjectDict):
    def __init__(self, transform):
        super(CODTransform, self).__init__()
        self.transform = transform

    def __getitem__(self, item):
        return self.post[self.transform(item)]

    def __setitem__(self, key, value):
        self.post[self.transform(key)] = value

    def _add_new(self, items, values):
        pass


class CODFunction(NestableObjectDict):
    def __init__(self, fnc, fnc_uses_lists=True):
        super(CODFunction, self).__init__()
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
#            raise KeyError('No cached values for %d items - %s' % (len(items), str(items)))

        if self.fnc_uses_lists:
            result = self._eval(items)
            return result
        else:
            return [self._eval(obj) for obj in items]

    def get_transformed_view(self, transform):
        def fnc(obj):
            return transform(self(obj))

        return fnc

class CODStore(NestableObjectDict):
    def __init__(self, name, dimensions, store, scope=None):
        super(CODStore, self).__init__()
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
        if isinstance(items, collections.Iterable):
            [dict.__setitem__(self, item, value) for item, value in zip(items, values)]
        else:
            dict.__setitem__(self, items, values)

        if self.max_save_buffer_size is not None and len(self) > self.max_save_buffer_size:
            self.sync()

    @property
    def storage(self):
        return self.store.storage

    def sync(self, flush_storable=True):
        storable = [ (key.idx[self.storage], value)
            for key, value in self.iteritems()
                if type(key) is not tuple and len(key.idx) > 0]

        if len(storable) > 0:
            storable_sorted = sorted(storable, key=lambda x: x[0])
            storable_keys = [x[0] for x in storable_sorted]
            storable_values = [x[1] for x in storable_sorted]
            self.store.set_list_value(self.scope, storable_keys, storable_values)

            if not flush_storable:
                non_storable = { key : value for key, value in self.iteritems() if len(key.idx) == 0 }
                self.clear()
                self.update(non_storable)
            else:
                self.clear()
        else:
            if flush_storable:
                self.clear()



    def flush_unstorable(self):
        storable = { key : value for key, value in self.iteritems() if len(key.idx) > 0 }
        self.clear()
        self.update(storable)

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
            return [values_sorted[idx] for idx in loadable_idxs]
        else:
            return []

    def _basetype(self, item):
        if type(item) is tuple:
            return item[0].content_class
        elif hasattr(item, 'base_cls'):
            return item.base_cls
        else:
            return type(item)

class CODMultiStore(CODStore):
    def __init__(self, store_name, name, dimensions, scope):
        super(CODStore, self).__init__()
        self.name = name
        self.dimensions = dimensions
        self.store_name = store_name

        if scope is None:
            self.scope = self
        else:
            self.scope = scope

        self.cod_stores = {}
        self.update_nod_stores()

    @property
    def storages(self):
        if hasattr(self.scope, 'idx'):
            return self.scope.idx.keys()
        else:
            return []

    def flush_unstorable(self):
        if len(self.storages) != len(self.cod_stores):
            self.update_nod_stores()

        if len(self.cod_stores) == 0:
            return None

        [ store.flush_unstorable() for store in self.cod_stores.values() ]


    def sync(self, flush_storable=True):
        if len(self.storages) != len(self.cod_stores):
            self.update_nod_stores()

        if len(self.cod_stores) == 0:
            return None

        [ store.sync(flush_storable) for store in self.cod_stores.values() ]

    def add_nod_store(self, storage):
        self.cod_stores[storage] = CODStore(self.name, self.dimensions, getattr(storage, self.store_name), self.scope)

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

        pass

    def _get(self, item):
        if len(self.storages) != len(self.cod_stores):
            self.update_nod_stores()

        if len(self.cod_stores) == 0:
            return None

        results = dict()
        for s in self.cod_stores:
            results[s] = self.cod_stores[s][item]

        for s, result in results.iteritems():
            if result is not None:
                return result

        return None


    def _get_list(self, items):
        if len(self.storages) != len(self.cod_stores):
            self.update_nod_stores()

        if len(self.cod_stores) == 0:
            return [None] * len(items)

        results_list = dict()
        for s in self.cod_stores:
            results_list[s] = self.cod_stores[s][items]

        output = [None] * len(items)
        for s, results in results_list.iteritems():
            output = [item if item is not None or results is None else result
                     for item, result in zip(output, results) ]

        return output

class CODUnwrapTuple(NestableObjectDict):
    def __init__(self):
        super(CODUnwrapTuple, self).__init__()

    def __getitem__(self, items):
        if isinstance(items, collections.Iterable):
            return self.post([value[0].load(value[1])
                if type(value) is tuple else value for value in items])
        else:
            if type(items) is tuple:
                items = items[0].load(items[1])

            return self.post[items]

    def __setitem__(self, key, value):
        self.post[key] = value




class CODWrap(NestableObjectDict):
    def __init__(self, post):
        super(CODWrap, self).__init__()
        self.post = post

    def __getitem__(self, items):
        return self.post[items]

    def __setitem__(self, key, value):
        self.post[key] = value

class CODBuffer(NestableObjectDict):
    """
    Implements a dict with a buffer that reads sequentially ahead
    """

class CODCache(NestableObjectDict):
    """
    Implements a dict with intelligent caching of limited size
    """