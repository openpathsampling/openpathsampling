from collections import OrderedDict
import weakref

__author__ = 'Jan-Hendrik Prinz'


class Cache(object):
    @property
    def count(self):
        return len(self), 0

    @property
    def size(self):
        return -1, -1

    def __str__(self):
        size = self.count
        maximum = self.size
        return '%s(%d/%d of %s/%s)' % (
            self.__class__.__name__,
            size[0], size[1],
            'Inf' if maximum[0] < 0 else str(maximum[0]),
            'Inf' if maximum[1] < 0 else str(maximum[1])
        )

    def __getitem__(self, item):
        raise KeyError("No items")

    def __setitem__(self, key, value):
        pass

    def get(self, item, default=None):
        try:
            return self[item]
        except KeyError:
            return default

    def transfer(self, old_cache):
        size = self.size
        if size[0] == -1 or size[1] == -1:
            for key in reversed(list(old_cache)):
                try:
                    self[key] = old_cache[key]
                except KeyError:
                    pass
        else:
            for key in reversed(list(old_cache)[0:size[0] + size[1]]):
                try:
                    self[key] = old_cache[key]
                except KeyError:
                    pass

        return self


class NoCache(Cache):
    def __init__(self):
        super(NoCache, self).__init__()

    def __getitem__(self, item):
        raise KeyError('No Cache has no items')

    def __contains__(self, item):
        return False

    def __setitem__(self, key, value):
        pass

    def __iter__(self):
        return iter([])

    @property
    def count(self):
        return 0, 0

    @property
    def size(self):
        return 0, 0

    def items(self):
        return []

    def transfer(self, old_cache):
        return self

    def clear(self):
        pass


class MaxCache(dict, Cache):
    def __init__(self):
        super(MaxCache, self).__init__()
        Cache.__init__(self)

    @property
    def count(self):
        return len(self), 0

    @property
    def size(self):
        return -1, 0


class LRUCache(Cache):
    """
    Implements a simple Least Recently Used Cache

    Very simple using collections.OrderedDict. The size can be change during
    run-time
    """

    def __init__(self, size_limit):
        super(LRUCache, self).__init__()
        self._size_limit = size_limit
        self._cache = OrderedDict()

    @property
    def count(self):
        return len(self._cache), 0

    @property
    def size(self):
        return self.size_limit, 0

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

    def __reversed__(self):
        return reversed(self._cache)

    def __getitem__(self, item):
        obj = self._cache.pop(item)
        self._cache[item] = obj
        return obj

    def __setitem__(self, key, value, **kwargs):
        self._cache[key] = value
        self._check_size_limit()

    def _check_size_limit(self):
        while len(self._cache) > self.size_limit:
            self._cache.popitem(last=False)

    def __contains__(self, item):
        return item in self._cache

    def clear(self):
        self._cache.clear()

    def __len__(self):
        return len(self._cache)


class WeakLRUCache(Cache):
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

        super(WeakLRUCache, self).__init__()
        self._size_limit = size_limit
        self.weak_type = weak_type

        if weak_type == 'value':
            self._weak_cache = weakref.WeakValueDictionary()
        elif weak_type == 'key':
            self._weak_cache = weakref.WeakKeyDictionary()
        else:
            raise ValueError("weak_type must be either 'key' or 'value'")

        self._cache = OrderedDict()

    @property
    def count(self):
        return len(self._cache), len(self._weak_cache)

    @property
    def size(self):
        return self._size_limit, -1

    def clear(self):
        self._cache.clear()
        self._weak_cache.clear()

    @property
    def size_limit(self):
        return self._size_limit

    def __getitem__(self, item):
        try:
            obj = self._cache.pop(item)
            self._cache[item] = obj
            return obj
        except KeyError:
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

    def get_silent(self, item):
        """
        Return item from the dict if it exists, None otherwise without reordering the LRU

        Parameters
        ----------
        item : object
            the item index to be retrieved from the cache

        Returns
        -------
        object of None
            the requested object if it exists else `None`
        """
        if item is None:
            return None

        try:
            return self._cache[item]
        except KeyError:
            try:
                return self._weak_cache[item]
            except KeyError:
                return None

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

    def __len__(self):
        return len(self._cache) + len(self._weak_cache)

    def __iter__(self):
        for key in self.keys():
            yield key

    def __reversed__(self):
        for key in reversed(self._weak_cache):
            yield key

        for key in reversed(self._cache):
            yield key


class WeakValueCache(weakref.WeakValueDictionary, Cache):
    """
    Implements a cache that keeps weak references to all elements
    """

    def __init__(self, *args, **kwargs):
        weakref.WeakValueDictionary.__init__(self, *args, **kwargs)
        Cache.__init__(self)

    @property
    def count(self):
        return 0, len(self)

    @property
    def size(self):
        return 0, -1


class WeakKeyCache(weakref.WeakKeyDictionary, Cache):
    """
    Implements a cache that keeps weak references to all elements
    """

    @property
    def count(self):
        return 0, len(self)

    @property
    def size(self):
        return 0, -1


class LRUChunkLoadingCache(Cache):
    """
    Implements a cache that keeps weak references to all elements

    In addition it uses a simple Least Recently Used Cache to make sure a portion
    of the last used elements are still present. Usually this number is 100.

    """

    def __init__(self, chunksize=100, max_chunks=100, variable=None):

        super(LRUChunkLoadingCache, self).__init__()
        self.max_chunks = max_chunks
        self.chunksize = chunksize
        self.variable = variable

        self._chunkdict = OrderedDict()
        self._firstchunk = 0
        self._lastchunk = []
        if variable is not None:
            self._size = len(self.variable)
        else:
            self._size = 0

        self._lastchunk_idx = self._size / self.chunksize

    @property
    def count(self):
        return sum(map(len, self._chunkdict)), 0

    @property
    def size(self):
        return self.max_chunks * self.chunksize, 0

    def clear(self):
        self._chunkdict.clear()

    def load_chunk(self, chunk_idx):
        if chunk_idx >= self._lastchunk_idx:
            self._size = len(self.variable)
            self._lastchunk_idx = self._size / self.chunksize

        if chunk_idx < self._lastchunk_idx:
            if chunk_idx not in self._chunkdict:
                # chunk not cached, load full
                left = chunk_idx * self.chunksize
                right = min(len(self.variable), left + self.chunksize)
                self._chunkdict[chunk_idx] = self.variable[left:right]

                self._check_size_limit()

            elif len(self._chunkdict[chunk_idx]) < self.chunksize:
                # incomplete chunk, load rest
                chunk = self._chunkdict[chunk_idx]

                left = chunk_idx * self.chunksize + len(chunk)
                right = min(len(self.variable), (chunk_idx + 1) * self.chunksize)

                if right > left:
                    chunk.append(self.variable[left:right])

    def _update_chunk_order(self, chunk_idx):
        if chunk_idx != self._firstchunk:
            chunk = self._chunkdict[chunk_idx]
            del self._chunkdict[chunk_idx]
            self._chunkdict[chunk_idx] = chunk
            self._firstchunk = chunk_idx

    def __getitem__(self, item):
        chunksize = self.chunksize
        chunk_idx = item / chunksize
        if chunk_idx in self._chunkdict:
            try:
                obj = self._chunkdict[chunk_idx][item % chunksize]
                self._update_chunk_order(chunk_idx)
                return obj
            except KeyError:
                pass

        self.load_chunk(chunk_idx)

        return self._chunkdict[chunk_idx][item % chunksize]

    def load_max(self):
        map(self.load_chunk,
            range(0, min(1 + (len(self.variable) - 1) / self.chunksize, self.max_chunks)))

    def __setitem__(self, key, value, **kwargs):
        chunk_idx = key / self.chunksize
        if chunk_idx in self._chunkdict:
            chunk = self._chunkdict[chunk_idx]
        else:
            chunk = []
            self._chunkdict[chunk_idx] = chunk

        left = chunk_idx * self.chunksize + len(chunk)
        right = key - 1

        if right > left:
            chunk.append(self.variable[left:right])

        chunk.append(value)

        self._update_chunk_order(chunk_idx)
        self._check_size_limit()

    def get_silent(self, item):
        """
        Return item from the dict if it exists, None otherwise without reordering the LRU

        Parameters
        ----------
        item : object
            the item index to be retrieved from the cache

        Returns
        -------
        object of None
            the requested object if it exists else `None`
        """
        return self.get(item)

    def _check_size_limit(self):
        if len(self._chunkdict) > self.max_chunks:
            self._chunkdict.popitem(last=False)

    def __contains__(self, item):
        return any(item in chunk for chunk in self._chunkdict)

    def keys(self):
        return sum(map(lambda x: x.keys(), self._chunkdict))

    def values(self):
        return sum(map(lambda x: x.values(), self._chunkdict))

    def __len__(self):
        return sum(map(len, self._chunkdict))

    def __iter__(self):
        for chunk in self._chunkdict.itervalues():
            for key in chunk.keys():
                yield key

    def __reversed__(self):
        for chunk in reversed(self._chunkdict.values()):
            for key in reversed(chunk.keys()):
                yield key
