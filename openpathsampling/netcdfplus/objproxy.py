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
    __slots__ = ['_subject', '_idx', '_store', '__weakref__']

    def __init__(self, store, idx):
        self._idx = idx
        self._store = store
        self._subject = None

    @property
    def __subject__(self):
        if self._subject is not None:
            obj = self._subject()
            if obj is not None:
                return obj

        ref = self._load_()

        if ref is None:
            return None

        self._subject = weakref.ref(ref)
        return ref

    def __eq__(self, other):
        if self is other:
            return True
        elif type(other) is LoaderProxy:
            if self._idx == other._idx and self._store is other._store:
                return True
        elif self.__subject__ is other:
            return True

        return False

    @property
    def __class__(self):
        return self._store.content_class

    def __getattr__(self, item):
        return getattr(self.__subject__, item)

    def _load_(self):
        """
        Call the loader and get the referenced object
        """
        return self._store[self._idx]
