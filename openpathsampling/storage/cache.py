__author__ = 'jan-hendrikprinz'

from collections import OrderedDict
import weakref

class LRUCache(OrderedDict):
    """
    Implements a simple Least Recently Used Cache

    Very simple using collections.OrderedDict. The size can be change during
    run-time
    """
    def __init__(self, size_limit):
        self._size_limit = size_limit
        OrderedDict.__init__(self)
        self._check_size_limit()

    @property
    def size_limit(self):
        return self._size_limit

    @size_limit.setter
    def size_limit(self, new_size):
        if new_size < self.size_limit:
          self._check_size_limit()

        self._size_limit = new_size

    def __setitem__(self, key, value, **kwargs) :
        OrderedDict.__setitem__(self, key, value)
        self._check_size_limit()

    def _check_size_limit(self):
        if self.size_limit is not None:
            while len(self) > self.size_limit:
                self.popitem(last=False)

class WeakLRUCache(OrderedDict):
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
        OrderedDict.__init__(self)

        self._size_limit = size_limit
        self._weak_cache = weakref.WeakValueDictionary()

    def __str__(self):
        return '%s(%d[%d])' % (self.__class__.__name__, len(self), len(self._weak_cache))

    @property
    def size_limit(self):
        return self._size_limit

    def __getitem__(self, item):
        try:
            return OrderedDict.__getitem__(self, item)
        except(KeyError):
            return self._weak_cache[item]

    @size_limit.setter
    def size_limit(self, new_size):
        if new_size < self.size_limit:
          self._check_size_limit()

        self._size_limit = new_size

    def __setitem__(self, key, value, **kwargs) :
        OrderedDict.__setitem__(self, key, value)
        self._check_size_limit()

    def _check_size_limit(self):
        if self.size_limit is not None:
            while len(self) > self.size_limit:
                self._weak_cache.__setitem__(*self.popitem(last=False))

    def __contains__(self, item):
        return OrderedDict.__contains__(self, item) or item in self._weak_cache

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