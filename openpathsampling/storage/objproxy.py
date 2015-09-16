'''

@author: JH Prinz
'''

from openpathsampling.base import StorableObject

import weakref

#=============================================================================
# SIMULATION CONFIGURATION
#=============================================================================

class LoaderProxy(StorableObject):
    """
    A proxy that loads an underlying object if attributes are accessed
    """

    def __init__(self):
        super(LoaderProxy, self).__init__()

    @property
    def __subject__(self):
        return self._load_()

    @property
    def __class__(self):
        store, idx = self.idx.iteritems().next()
        return store.content_class

    def __getattr__(self, item):
        return getattr(self.__subject__, item)

    def _load_(self):
        """
        Call the loader and get the referenced object
        """
        store, idx = self.idx.iteritems().next()
        obj = store[idx]

        return obj

        # print 'loaded %s[%d] : %s' % (store.content_class.__name__, idx, self.__subject__)

class DelayedLoaderProxy(StorableObject):
    """
    A proxy that loads an underlying object if attributes are accessed
    """

    def __init__(self):
        super(DelayedLoaderProxy, self).__init__()
        self._subject_ref_ = None

    @property
    def __subject__(self):
        if self._subject_ref_ is None:
            self._subject_ref_ = self._load_()

        return self._subject_ref_

    @property
    def __class__(self):
        return self.store.content_class

    def __getattr__(self, item):
        return getattr(self.__subject__, item)

    def _load_(self):
        """
        Call the loader and get the referenced object
        """
        store, idx = self.idx.iteritems().next()
        return store.get(idx) # .load would just get another Proxy

        # print 'loaded %s[%d] : %s' % (store.content_class.__name__, idx, self.__subject__)

    def unload(self):
        """
        Unload the referenced object to free memory
        """
        self._subject_ref_ = None

class WeakLoaderProxy(DelayedLoaderProxy):
    """
    A proxy that loads an underlying object with weak reference if attributes are accessed
    """

    def __init__(self):
        super(WeakLoaderProxy, self).__init__()
        self._subject_weak_ = None

    @property
    def __subject__(self):
        if self._subject_ref_ is None:
            if self._subject_weak_ is None:
                obj = self.load()
                self._subject_weak_ = weakref.ref(obj)
            else:
                obj = self._subject_weak_()
                if obj is None:
                    obj = self._load_()
                    self._subject_weak_ = weakref.ref(obj)

            return obj
        else:
            return self._subject_ref_

    def stick(self):
        self._subject_ref_ = self.__subject__

    def unstick(self):
        self._subject_weak_ = None

    def unload(self):
        self._subject_ref_ = None
        self._subject_weak_ = None