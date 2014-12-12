from opentis.storage import ObjectStorage
from wrapper import savecache, loadcache
from opentis.sample import Sample


class SampleStorage(ObjectStorage):
    def __init__(self, storage):
        super(SampleStorage, self).__init__(storage, Sample)

    @savecache
    def save(self, sample, idx=None):
        """
        Add the current state of the sample in the database. If nothing has changed then the sample gets stored using the same snapshots as before. Saving lots of diskspace

        Parameters
        ----------
        sample : Sample()
            the sample to be saved
        idx : int or None
            if idx is not None the index will be used for saving in the storage. This might overwrite already existing trajectories!

        Notes
        -----
        This also saves all contained frames in the sample if not done yet.
        A single Sample object can only be saved once!
        """

        if idx is not None:
            self.storage.trajectory.save(sample.trajectory)
            self.set_object('sample_trajectory', idx, sample.trajectory)

            self.storage.ensemble.save(sample.ensemble)
            self.set_object('sample_ensemble', idx, sample.ensemble)

            self.storage.pathmover.save(sample.mover)
            self.set_object('sample_mover', idx, sample.mover)

            self.storage.movedetails.save(sample.details)
            self.set_object('sample_details', idx, sample.details)

            self.save_variable('sample_step', idx, sample.time)

    @loadcache
    def load(self, idx, momentum = True):
        '''
        Return a sample from the storage

        Parameters
        ----------
        idx : int
            index of the sample (counts from 1)

        Returns
        -------
        sample : Sample
            the sample
        '''
        trajectory_idx = int(self.storage.variables['sample_trajectory_idx'][idx])
        ensemble_idx = int(self.storage.variables['sample_ensemble_idx'][idx])
        mover_idx = int(self.storage.variables['sample_mover_idx'][idx])
        details_idx = int(self.storage.variables['sample_details_idx'][idx])
        step=self.load_variable('sample_step', idx)


        obj = Sample(
            trajectory=self.storage.trajectory.load(trajectory_idx, lazy=True),
            mover=self.storage.pathmover.load(mover_idx, lazy=True),
            ensemble=self.storage.ensemble.load(ensemble_idx),
            details=self.storage.movedetails.load(details_idx),
            step=step
        )

        obj.idx[self.storage] = idx

        return obj

    def by_ensemble(self, ensemble):
        return [ sample for sample in self.iterator() if sample.ensemble == ensemble ]

    def _init(self):
        """
        Initialize the associated storage to allow for sample storage

        """
        super(SampleStorage, self)._init()

        # New short-hand definition
        self.init_variable('sample_trajectory_idx', 'index')
        self.init_variable('sample_ensemble_idx', 'index')
        self.init_variable('sample_mover_idx', 'index')
        self.init_variable('sample_details_idx', 'index')
        self.init_variable('sample_step', 'index')