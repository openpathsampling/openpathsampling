__author__ = 'jan-hendrikprinz'

from collections import OrderedDict
import weakref

class LRUCache(object):
    """
    Implements a simple Least Recently Used Cache

    Very simple using collections.OrderedDict. The size can be change during
    run-time
    """
    def __init__(self, size_limit):
        self._size_limit = size_limit
        self._cache = OrderedDict()

    @property
    def size_limit(self):
        return self._size_limit

    @size_limit.setter
    def size_limit(self, new_size):
        if new_size < self.size_limit:
          self._check_size_limit()

        self._size_limit = new_size

    def __iter__(self):
        return iter(self._cache)

    def __getitem__(self, item):
        obj = self._cache.pop(item)
        self._cache[item] = obj
        return obj

    def __setitem__(self, key, value, **kwargs) :
        self._cache[key] = value
        self._check_size_limit()

    def _check_size_limit(self):
        while len(self._cache) > self.size_limit:
            self._cache.popitem(last=False)

    def __str__(self):
        return '%s(%d)' % (self.__class__.__name__, len(self._cache))

    def __contains__(self, item):
        return item in self._cache


class WeakLRUCache(object):
    """
    Implements a cache that keeps weak references to all elements

    In addition it uses a simple Least Recently Used Cache to make sure a portion
    of the last used elements are still present. Usually this number is 100.

    """
    def __init__(self, size_limit=100, weak_type='value'):
        """
        Parameters
        ----------
        size_limit : int
            integer that defines the size of the LRU cache. Default is 100.
        """
        self._size_limit = size_limit
        self.weak_type = weak_type

        if weak_type == 'value':
            self._weak_cache = weakref.WeakValueDictionary()
        elif weak_type == 'key':
            self._weak_cache = weakref.WeakKeyDictionary()

        self._cache = OrderedDict()

    def __str__(self):
        return '%s(%d[%d])' % (self.__class__.__name__, len(self._cache), len(self._weak_cache))

    @property
    def size_limit(self):
        return self._size_limit

    def __getitem__(self, item):
        try:
            obj = self._cache.pop(item)
            self._cache[item] = obj
            return obj
        except(KeyError):
            obj = self._weak_cache[item]
            del self._weak_cache[item]
            self._cache[item] = obj
            self._check_size_limit()
            return obj

    @size_limit.setter
    def size_limit(self, new_size):
        if new_size < self.size_limit:
          self._check_size_limit()

        self._size_limit = new_size

    def __setitem__(self, key, value, **kwargs):
        try:
            self._cache.pop(key)
        except KeyError:
            pass

        self._cache[key] = value
        self._check_size_limit()


    def _check_size_limit(self):
        if self.size_limit is not None:
            while len(self._cache) > self.size_limit:
                self._weak_cache.__setitem__(*self._cache.popitem(last=False))

    def __contains__(self, item):
        return item in self._cache or item in self._weak_cache

    def keys(self):
        return self._cache.keys() + self._weak_cache.keys()

    def values(self):
        return self._cache.values() + self._weak_cache.values()

class WeakCache(weakref.WeakValueDictionary):
    """
    Implements a cache that keeps weak references to all elements
    """

class WeakLimitCache(dict):
    """
    Implements a cache that keeps weak references to all elements

    In addition it uses a simple Least Recently Used Cache to make sure a portion
    of the last used elements are still present. Usually this number is 100.

    """
    def __init__(self, size_limit=100):
        """
        Parameters
        ----------
        size_limit : int
            integer that defines the size of the LRU cache. Default is 100.
        """
        dict.__init__(self)

        self._size_limit = size_limit
        self._weak_cache = weakref.WeakValueDictionary()

    @property
    def size_limit(self):
        return self._size_limit

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except(KeyError):
            return self._weak_cache[item]

    @size_limit.setter
    def size_limit(self, new_size):
        if new_size < self.size_limit:
          self._check_size_limit()

        self._size_limit = new_size

    def __setitem__(self, key, value, **kwargs) :
        if self.size_limit is not None:
            if len(self) == self.size_limit:
                self._weak_cache[key] = value
            else:
                dict.__setitem__(self, key, value)

    def _check_size_limit(self):
        if self.size_limit is not None:
            while len(self) > self.size_limit:
                self._weak_cache.__setitem__(*self.popitem())