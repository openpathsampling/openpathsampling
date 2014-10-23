from object_storage import ObjectStorage, addcache_save, addcache

class Sample(object):
    """
    A Move is the return object from a PathMover and contains all information about the move, initial trajectories,
    new trajectories (both as references). Might move several trajectories at a time (swapping)

    Notes
    -----
    Should contain inputs/outputs and success (accepted/rejected) as well as probability to succeed.
    """

    cls = 'path'

    def __init__(self, trajectory=None,  mover=None, ensemble=None, details=None):
        self.idx = dict()

        self.mover = mover
        self.ensemble = ensemble
        self.trajectory = trajectory
        self.details = details

    def __call__(self):
        return self.trajectory

class SampleStorage(ObjectStorage):
    def __init__(self, storage):
        super(SampleStorage, self).__init__(storage, Sample)

    @addcache_save
    def save(self, origin, idx=None):
        """
        Add the current state of the origin in the database. If nothing has changed then the origin gets stored using the same snapshots as before. Saving lots of diskspace

        Parameters
        ----------
        origin : Sample()
            the origin to be saved
        idx : int or None
            if idx is not None the index will be used for saving in the storage. This might overwrite already existing trajectories!

        Notes
        -----
        This also saves all contained frames in the origin if not done yet.
        A single Sample object can only be saved once!
        """

        if idx is not None:
            storage = self.storage

            self.storage.trajectory.save(origin.trajectory)
            self.save_object('origin_trajectory', idx, origin.trajectory)

            self.storage.ensemble.save(origin.ensemble)
            self.save_object('origin_ensemble', idx, origin.ensemble)

            self.storage.pathmover.save(origin.mover)
            self.save_object('origin_mover', idx, origin.mover)

            self.storage.movedetails.save(origin.details)
            self.save_object('origin_details', idx, origin.details)

    @addcache
    def load(self, idx, momentum = True):
        '''
        Return a origin from the storage

        Parameters
        ----------
        idx : int
            index of the origin (counts from 1)

        Returns
        -------
        origin : Sample
            the origin
        '''
        trajectory_idx = self.storage.variables['origin_trajectory_idx'][idx]
        ensemble_idx = self.storage.variables['origin_ensemble_idx'][idx]
        mover_idx = self.storage.variables['origin_mover_idx'][idx]
        details_idx = self.storage.variables['origin_details_idx'][idx]

        obj = Sample(
            trajectory=self.storage.trajectory.load(trajectory_idx, lazy=True),
            mover=self.storage.pathmover.load(mover_idx, lazy=True),
            ensemble=self.storage.ensemble.load(ensemble_idx),
            details=self.storage.movedetails.load(details_idx)
        )

        return obj

    def _init(self):
        """
        Initialize the associated storage to allow for origin storage

        """
        super(SampleStorage, self)._init()

        # New short-hand definition
        self.init_variable('origin_trajectory_idx', 'u4')
        self.init_variable('origin_ensemble_idx', 'u4')
        self.init_variable('origin_mover_idx', 'u4')
        self.init_variable('origin_details_idx', 'u4')
