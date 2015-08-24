'''

@author: JH Prinz
'''

from openpathsampling.base import StorableObject

#=============================================================================
# SIMULATION CONFIGURATION
#=============================================================================

class DelayedLoaderProxy(StorableObject):
    """
    A proxy that loads an underlying object if attributes are accessed
    """

    # Class variables to store the global storage and the system context
    # describing the system to be saved as snapshots
    # Hopefully these class member variables will not be needed any longer
    engine = None

    def __init__(self):
        super(DelayedLoaderProxy, self).__init__()
        self.__subject__ = None

    @property
    def __class__(self):
        return self.store.content_class

    def __getattr__(self, item):
        if self.__subject__ is None:
            self.load()

        return getattr(self.__subject__, item)

    def load(self):
        """
        Call the loader and get the referenced object
        """
        if self.__subject__ is None:
            store, idx = self.idx.iteritems().next()
            self.__subject__ = store.get(idx) # .load would just get another Proxy

            # print 'loaded %s[%d] : %s' % (store.content_class.__name__, idx, self.__subject__)

    def unload(self):
        """
        Unload the referenced object to free memory
        """
        self.__subject__ = None