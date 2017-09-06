import collections
import numpy as np

from .proxy import LoaderProxy

__author__ = 'Jan-Hendrik Prinz'


class ChainDict(object):
    """
    A dict-like structure with a logic to fill missing values from other dicts

    The default ChainDict requires a list of keys. If you want to allow also
    single values you need to add a ChainDict that interprets single and
    iterables like ExpandSingle.

    The default for unknown keys is None. This is necessary to pass on what
    is unknown. Everything that is not known in the current instance is
    passed on to .post

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

    First time the function is called and the cache filled. Second time
    the cache is used.

    >>> combo_cd[[1,2]]
    [1,4]

    Attributes
    ----------
    _post : ChainDict
        the ChainDict to be called when this instance cannot evaluate given keys

    """

    def __init__(self):
        self._post = None
        self._iterables = (list, tuple, set, frozenset)
        self._singles = ()

    def __getitem__(self, items):
        # first apply the own _get functions to compute
        # print 'get %d items using %s' % (len(items), self.__class__.__name__)

        results = self._get_list(items)

        if self._post is not None:
            nones = [obj[0] for obj in zip(items, results) if obj[1] is None]
            if len(nones) == 0:
                return results
            else:
                rep = list(self._post[nones])
                self._set_list(nones, rep)

                it = iter(rep)
                return [next(it) if p[1] is None else p[1]
                        for p in zip(items, results)]

        return results

    __call__ = __getitem__

    def __setitem__(self, key, value):
        self._set_list(key, value)

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
        [self._set(item, value) for item, value in zip(items, values)
         if values is not None]

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
        return ' > '.join(
            map(lambda x: x.__class__.__name__, self.passing_chain))


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

    __call__ = __getitem__

    def __setitem__(self, key, value):
        self._post[key] = value


class MergeNumpy(ChainDict):
    """Wrap returned values in a numpy array
    """

    def __getitem__(self, items):
        return np.array(self._post[items])

    __call__ = __getitem__

    def __setitem__(self, key, value):
        self._post[key] = value


class ExpandSingle(ChainDict):
    """
    Iterables will be unrolled and passed as a list
    """
    def __init__(self, key_class):
        super(ExpandSingle, self).__init__()
        self.key_class = key_class

    def __getitem__(self, items):
        if isinstance(items, self.key_class):
            return self._post[[items]][0]
        # elif type(items) is LoaderProxy:
        #     return self._post[[items]][0]
        elif hasattr(items.__class__, '__iter__'):
            try:
                _ = len(items)
            except AttributeError:
                # no length means unbound iterator and we cannot handle these
                raise AttributeError(
                    'Iterators that do not have __len__ implemented are not '
                    'supported. You can wrap your iterator in list() if you '
                    'know that it will finish.')
            try:
                return self._post[items.as_proxies()]
            except AttributeError:
                # turn possible iterators into a list if possible
                return self._post[list(items)]

        else:
            return self._post[[items]][0]

    __call__ = __getitem__

    def __setitem__(self, items, values):
        if isinstance(items, self.key_class):
            self._post[[items]] = [values]
        elif type(items) is LoaderProxy:
            self._post[[items]] = [values]
        elif hasattr(items, '__iter__'):
            try:
                _ = len(items)
            except AttributeError:
                # no length means unbound iterator and we cannot handle these
                raise AttributeError(
                    'Iterators that do not have __len__ implemented are not '
                    'supported. You can wrap your iterator in list() if you '
                    'know that it will finish.')
            try:
                self._post[items.as_proxies()] = values
            except AttributeError:
                # turn possible iterators into a list if possible
                self._post[list(items)] = values

        else:
            self._post[[items]] = values


class Function(ChainDict):
    """
    Uses a regular function to evaluate given keys.

    This works effective like a function called with square brackets
    """
    def __init__(
            self,
            fnc,
            requires_lists=True,
            scalarize_numpy_singletons=False):
        """
        Parameters
        ----------
        fnc : function
            the function to be evaluated to return values to keys
        requires_lists : bool
            if `True` we assume that it is faster to pass lists to this
            function instead of evaluating each key separately
        scalarize_numpy_singletons : bool
            if `True` eventual numpy objects that have length one in their
            last dimension, will be flattened by the last dimension. This
            is often useful if you have function that will by default return
            a list of results. In case your function does so, you can
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
                results = map(lambda x: x.reshape(x.shape[:-1]), results)

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
    Return Values from a cache filled from underlying CDs
    """
    def __init__(self, cache):
        """
        Parameters
        ----------
        cache : :obj:`openpathsampling.netcdfplus.cache.Cache` or dict
            the cache to be used to store the data
        """
        super(CacheChainDict, self).__init__()
        self.cache = cache

    def _contains(self, item):
        return item in self.cache

    def _get(self, item):
        if item is None:
            return None

        try:
            return self.cache[item]
        except KeyError:
            return None

    def _set(self, item, value):
        self.cache[item] = value


class ReversibleCacheChainDict(CacheChainDict):
    """
    Return Values from a cache filled from the underlying CD
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

    def _get(self, item):
        if item is None:
            return None

        try:
            return self.cache[item]
        except KeyError:
            if type(item) is not LoaderProxy:
                if self.reversible and item._reversed is not None:
                    try:
                        return self.cache[item._reversed]
                    except KeyError:
                        return None

            return None


class StoredDict(ChainDict):
    """
    ChainDict that has a store attached and returns existing store values
    """
    def __init__(self, value_store):
        """
        Parameters
        ----------
        value_store : storage.Variable
            the store that references the store variable to store the values
            by index
        """
        super(StoredDict, self).__init__()
        self.value_store = value_store

    def _set(self, item, value):
        return

    def _set_list(self, items, values):
        return

    def _get(self, item):
        return self.value_store.get(item)

    def _get_list(self, items):
        return list(map(self._get, items))

    def sync(self):
        pass

    def cache_all(self):
        pass

