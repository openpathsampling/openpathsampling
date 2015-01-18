from object_storage import ObjectStorage
from opentis.sample import Sample, SampleSet

class SampleStorage(ObjectStorage):
    def __init__(self, storage):
        super(SampleStorage, self).__init__(storage, Sample)

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

            self.save_variable('sample_replica', idx, sample.replica)

            self.storage.movedetails.save(sample.details)
            self.set_object('sample_details', idx, sample.details)

            self.save_variable('sample_step', idx, sample.time)

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
        replica_idx = int(self.storage.variables['sample_replica'][idx])
        details_idx = int(self.storage.variables['sample_details_idx'][idx])
        step=self.load_variable('sample_step', idx)


        obj = Sample(
            trajectory=self.storage.trajectory.load(trajectory_idx, lazy=True),
            replica=replica_idx,
            ensemble=self.storage.ensemble.load(ensemble_idx),
            details=self.storage.movedetails.load(details_idx),
            step=step
        )

        return obj

    def by_ensemble(self, ensemble):
        return [ sample for sample in self.iterator() if sample.ensemble == ensemble ]

    def _init(self):
        """
        Initialize the associated storage to allow for sample storage

        """
        super(SampleStorage, self)._init()

        # New short-hand definition
        self.init_variable('sample_trajectory_idx', 'index', chunksizes=(1, ))
        self.init_variable('sample_ensemble_idx', 'index', chunksizes=(1, ))
        self.init_variable('sample_replica', 'index', chunksizes=(1, ))
        self.init_variable('sample_details_idx', 'index', chunksizes=(1, ))
        self.init_variable('sample_step', 'index', chunksizes=(1, ))

class SampleSetStorage(ObjectStorage):

    def __init__(self, storage):
        super(SampleSetStorage, self).__init__(storage, SampleSet)

    def save(self, sampleset, idx=None):
        """
        Add the current state of the sampleset in the database. If nothing has changed then the sampleset gets stored using the same samples as before. Saving lots of diskspace

        Parameters
        ----------
        sampleset : Trajectory()
            the sampleset to be saved
        idx : int or None
            if idx is not None the index will be used for saving in the storage. This might overwrite already existing trajectories!

        Notes
        -----
        This also saves all contained frames in the sampleset if not done yet.
        A single Trajectory object can only be saved once!
        """

        # Check if all samples are saved
        map(self.storage.sample.save, sampleset)

        values = self.list_to_numpy(sampleset, 'sample')
        self.storage.variables['sampleset_sample_idx'][idx] = values


    def sample_indices(self, idx):
        '''
        Load sample indices for sampleset with ID 'idx' from the storage

        ARGUMENTS

        idx (int) - ID of the sampleset

        Returns
        -------
        list of int - sampleset indices
        '''

        # get the values
        values = self.storage.variables['sampleset_sample_idx'][idx]

        # typecast to integer
        return self.list_from_numpy(values, 'index')

    def load(self, idx, lazy = None):
        '''
        Return a sampleset from the storage

        Parameters
        ----------
        idx : int
            index of the sampleset (counts from 0)

        Returns
        -------
        sampleset : Trajectory
            the sampleset
        '''

        values = self.storage.variables['sampleset_sample_idx'][idx]

        # typecast to sample
        samples = self.list_from_numpy(values, 'sample')
        sampleset = SampleSet(samples)

        return sampleset

    def _init(self, units=None):
        """
        Initialize the associated storage to allow for sampleset storage

        """
        super(SampleSetStorage, self)._init()

        self.init_variable('sampleset_sample_idx', 'index', 'sampleset',
            description="sampleset[sampleset][frame] is the sample index (0..nspanshots-1) of frame 'frame' of sampleset 'sampleset'.",
            variable_length = True,
            chunksizes=(1024, )
        )