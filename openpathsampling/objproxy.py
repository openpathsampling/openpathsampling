'''

@author: JH Prinz
'''

from openpathsampling.todict import OPSObject

#=============================================================================
# SIMULATION CONFIGURATION
#=============================================================================

class DelayedLoaderProxy(OPSObject):
    """
    Simulation snapshot. Contains references to a configuration and momentum
    """

    # Class variables to store the global storage and the system context
    # describing the system to be saved as snapshots
    # Hopefully these class member variables will not be needed any longer
    engine = None

    def __init__(self):
        """

        """
        super(DelayedLoaderProxy, self).__init__()
        self.__subject__ = None

    @property
    def __class__(self):
        return self.store.content_class

    def __getattr__(self, item):
        if self.__subject__ is None:
            store, idx = self.idx.iteritems().next()
            self.__subject__ = store.get(idx) # .load would just get another Proxy

#            print 'loaded %s[%d] : %s' % (store.content_class.__name__, idx, self.__subject__)

        return getattr(self.__subject__, item)