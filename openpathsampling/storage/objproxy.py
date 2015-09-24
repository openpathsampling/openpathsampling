'''

@author: JH Prinz
'''

from openpathsampling.base import StorableObject

import weakref

#=============================================================================
# SIMULATION CONFIGURATION
#=============================================================================

class LoaderProxy(object):
    """
    A proxy that loads an underlying object if attributes are accessed
    """


    def __init__(self, idx=None):
        if idx is None:
            idx = dict()
        super(LoaderProxy, self).__init__()
        self._idx = idx

    @property
    def __subject__(self):
        return self._load_()

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