from snapshot import Snapshot, Configuration, Momentum

from storage_utils import setstorage

from functools import wraps

class SnapshotStorage(object):

    def __init__(self, storage = None):
        self.storage = storage

    @wraps(setstorage)
    def save(self, snapshot, idx_configuration = None, idx_momentum = None):
        """
        Save positions, velocities, boxvectors and energies of current iteration to NetCDF file.

        Notes
        -----
        We need to allow for reversed snapshots to save memory. Would be nice
        """

        self.storage.configuration.save(snapshot.configuration, idx_configuration)
        self.storage.momentum.save(snapshot.momentum, idx_momentum)

    @wraps(setstorage)
    def load(self, idx_configuration = None, idx_momentum = None, reversed = False):
        '''
        Load a snapshot from the storage

        Parameters
        ----------
        idx : int
            index of the snapshot in the database 'idx' > 0

        Returns
        -------
        snapshot : Snapshot
            the snapshot
        '''

        #TODO: Check, for some reason some idx are given as numpy.in32 and netcdf4 is not compatible with indices given in this format!!!!!

        snapshot = Snapshot()
        snapshot.reversed = bool(reversed)

        if idx_configuration is not None:
            idx_c = int(idx_configuration)
            snapshot.configuration = self.storage.configuration.load(idx_c)
            print snapshot.configuration.idx

        if idx_momentum is not None:
            idx_m = int(idx_momentum)
            snapshot.momentum = self.storage.momentum.load(idx_m)

        return snapshot

    def _init(self):
        """
        Initialize the associated storage to allow for snapshot storage

        """