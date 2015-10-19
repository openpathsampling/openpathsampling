"""

@author: JH Prinz
"""

import weakref

# =============================================================================
# Loader Proxy
# =============================================================================

class LoaderProxy(object):
    """
    A proxy that loads an underlying object if attributes are accessed
    """

    def __init__(self, idx=None):
        if idx is None:
            idx = dict()
        super(LoaderProxy, self).__init__()
        self._idx = idx
        self._subject = None

    @property
    def __subject__(self):
        if self._subject is not None:
            obj = self._subject()
            if obj is not None:
                return obj

        ref = self._load_()
        self._subject = weakref.ref(ref)
        return ref

    def __eq__(self, other):
        if self is other:
            return True
        elif self.__subject__ is other:
            return True
        elif type(other) is LoaderProxy:
            store1, idx1 = self._idx.iteritems().next()
            store2, idx2 = other._idx.iteritems().next()

            if idx1 == idx2 and store1 is store2:
                return True

        return False

    @property
    def __class__(self):
        store, idx = self._idx.iteritems().next()
        return store.content_class

    def __getattr__(self, item):
        return getattr(self.__subject__, item)

    def _load_(self):
        """
        Call the loader and get the referenced object
        """
        store, idx = self._idx.iteritems().next()
        obj = store[idx]

        return obj
