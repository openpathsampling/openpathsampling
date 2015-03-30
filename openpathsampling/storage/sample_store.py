from object_storage import ObjectStore
from openpathsampling.sample import SampleSet, Sample

import time

class SampleStore(ObjectStore):
    def __init__(self, storage):
        super(SampleStore, self).__init__(storage, Sample, json=False, load_partial=True)

        self.set_variable_partial_loading('details', self.update_details)

    def save(self, sample, idx=None):
        if idx is not None:
            self.storage.trajectories.save(sample.trajectory)
            self.set_object('sample_trajectory', idx, sample.trajectory)

            self.storage.ensembles.save(sample.ensemble)
            self.set_object('sample_ensemble', idx, sample.ensemble)

            self.save_variable('sample_replica', idx, sample.replica)

            self.storage.movedetails.save(sample.details)
            self.set_object('sample_details', idx, sample.details)

            if sample.step is None:
                self.save_variable('sample_step', idx, -1)
            else:
                self.save_variable('sample_step', idx, sample.step)

    def load_empty(self, idx):
        trajectory_idx = int(self.storage.variables['sample_trajectory_idx'][idx])
        ensemble_idx = int(self.storage.variables['sample_ensemble_idx'][idx])
        replica_idx = int(self.storage.variables['sample_replica'][idx])
#        details_idx = int(self.storage.variables['sample_details_idx'][idx])
        step=self.load_variable('sample_step', idx)

        if step < 0:
            step = None


        obj = Sample(
            trajectory=self.storage.trajectories.load(trajectory_idx),
            replica=replica_idx,
            ensemble=self.storage.ensembles.load(ensemble_idx),
            step=step
        )

        del obj.details

        return obj

    @staticmethod
    def update_details(obj):
        storage = obj._origin

        idx = obj.idx[storage]
        details_idx = int(storage.variables['sample_details_idx'][idx])
        details = storage.movedetails.load(details_idx)

        obj.details = details

    def load(self, idx):
        '''
        Return a sample from the storage

        Parameters
        ----------
        idx : int
            index of the sample

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
            trajectory=self.storage.trajectories.load(trajectory_idx),
            replica=replica_idx,
            ensemble=self.storage.ensembles.load(ensemble_idx),
            details=self.storage.movedetails.load(details_idx),
            step=step
        )

        return obj

    def by_ensemble(self, ensemble):
        return [ sample for sample in self.iterator() if sample.ensemble == ensemble ]

    def _init(self):
        super(SampleStore, self)._init()

        # New short-hand definition
        self.init_variable('sample_trajectory_idx', 'index', chunksizes=(1, ))
        self.init_variable('sample_ensemble_idx', 'index', chunksizes=(1, ))
        self.init_variable('sample_replica', 'index', chunksizes=(1, ))
        self.init_variable('sample_details_idx', 'index', chunksizes=(1, ))
        self.init_variable('sample_step', 'index', chunksizes=(1, ))

class SampleSetStore(ObjectStore):

    def __init__(self, storage):
        super(SampleSetStore, self).__init__(storage, SampleSet, json=False)

    def save(self, sample_set, idx=None):
        # Check if all samples are saved
        map(self.storage.samples.save, sample_set)

        values = self.list_to_numpy(sample_set, 'sample')
        self.storage.variables['sampleset_sample_idx'][idx] = values

        self.storage.movepaths.save(sample_set.movepath)
        self.set_object('sampleset_movepath', idx, sample_set.movepath)


    def sample_indices(self, idx):
        '''
        Load sample indices for sample_set with ID 'idx' from the storage

        Parameters
        ----------
        idx : int
            ID of the sample_set

        Returns
        -------
        list of int
            list of sample indices
        '''

        # get the values
        values = self.storage.variables['sampleset_sample_idx'][idx]

        # typecast to integer
        return self.list_from_numpy(values, 'index')

    def load(self, idx):
        '''
        Return a sample_set from the storage

        Parameters
        ----------
        idx : int
            index of the sample_set (counts from 0)

        Returns
        -------
        sample_set
            the sample_set
        '''

        movepath_idx = int(self.storage.variables['sampleset_movepath_idx'][idx])

        values = self.storage.variables['sampleset_sample_idx'][idx]

        # typecast to sample
        samples = self.list_from_numpy(values, 'sample')
        sample_set = SampleSet(samples, movepath=self.storage.movepaths.load(movepath_idx))

        return sample_set

    def _init(self, units=None):
        """
        Initialize the associated storage to allow for sampleset storage

        """
        super(SampleSetStore, self)._init()

        self.init_variable('sampleset_sample_idx', 'index', 'sampleset',
            description="sampleset[sampleset][frame] is the sample index (0..nspanshots-1) of frame 'frame' of sampleset 'sampleset'.",
            variable_length = True,
            chunksizes=(1024, )
        )

        self.init_variable('sampleset_movepath_idx', 'index', chunksizes=(1, ))
