from object_storage import ObjectStore
from openpathsampling.sample import SampleSet, Sample
from openpathsampling.storage.object_storage import func_update_object

import time

class SampleStore(ObjectStore):
    def __init__(self, storage):
        super(SampleStore, self).__init__(storage, Sample, json=False, load_partial=True)

        self.set_variable_partial_loading('details', self.update_details)
        self.set_variable_partial_loading('parent', self.update_parent)
        self.set_variable_partial_loading('mover', self.update_mover)
        self.set_variable_partial_loading('ensemble', self.update_ensemble)

        self._cached_all = False


    def load_empty(self, idx):
        trajectory_idx = int(self.storage.variables['sample_trajectory_idx'][idx])
        replica_idx = int(self.storage.variables['sample_replica'][idx])
        bias=float(self.load_variable('sample_bias', idx))

        obj = Sample(
            trajectory=self.storage.trajectories[trajectory_idx],
            replica=replica_idx,
            bias=bias,
        )

        del obj.details
        del obj.ensemble
        del obj.mover
        del obj.parent

        return obj

    update_details = func_update_object('details', 'sample', 'details', '_details')
    update_parent = func_update_object('parent', 'sample', 'parent', 'samples')
    update_mover = func_update_object('mover', 'sample', 'pathmover', 'pathmovers')
    update_ensemble = func_update_object('ensemble', 'sample', 'ensemble', 'ensembles')

    def save(self, sample, idx=None):
        if idx is not None:
            self.storage.trajectories.save(sample.trajectory)
            self.save_object('sample_trajectory', idx, sample.trajectory)

            self.storage.ensembles.save(sample.ensemble)
            self.save_object('sample_ensemble', idx, sample.ensemble)

            self.save_variable('sample_replica', idx, sample.replica)
            self.save_object('sample_parent', idx, sample.parent)
            self.save_object('sample_details', idx, sample.details)
            self.save_variable('sample_bias', idx, sample.bias)
            self.save_object('sample_pathmover', idx, sample.mover)

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
        parent_idx = int(self.storage.variables['sample_parent_idx'][idx])
        details_idx = int(self.storage.variables['sample_details_idx'][idx])
        pathmover_idx = int(self.storage.variables['sample_pathmover_idx'][idx])
        bias=float(self.load_variable('sample_bias', idx))


        obj = Sample(
            trajectory=self.storage.trajectories[trajectory_idx],
            replica=replica_idx,
            ensemble=self.storage.ensembles[ensemble_idx],
            parent=self.storage.samples[parent_idx],
            details=self.storage._details[details_idx],
            bias=bias,
            mover=self.storage.pathmovers[pathmover_idx]
        )

        return obj

    def by_ensemble(self, ensemble):
        return [ sample for sample in self.iterator() if sample.ensemble == ensemble ]

    def _init(self, units=None):
        super(SampleStore, self)._init(units)

        # New short-hand definition
        self.init_variable('sample_trajectory_idx', 'index', chunksizes=(1, ))
        self.init_variable('sample_ensemble_idx', 'index', chunksizes=(1, ))
        self.init_variable('sample_replica', 'index', chunksizes=(1, ))
        self.init_variable('sample_parent_idx', 'index', chunksizes=(1, ))
        self.init_variable('sample_details_idx', 'index', chunksizes=(1, ))
        self.init_variable('sample_bias', 'float', chunksizes=(1, ))
        self.init_variable('sample_pathmover_idx', 'index', chunksizes=(1, ))

    def all(self):
        self.cache_all()
        return self

    def cache_all(self):
        """Load all samples as fast as possible into the cache

        """
        if not self._cached_all:
            idxs = range(len(self))
            trajectory_idxs = self.storage.variables['sample_trajectory_idx'][:]
            replica_idxs = self.storage.variables['sample_replica'][:]
            biass = self.storage.variables['sample_bias'][:]

            [ self.add_empty_to_cache(i,t,r,a) for i,t,r,a in zip(
                idxs,
                trajectory_idxs,
                replica_idxs,
                biass) ]

            self._cached_all = True


    def add_empty_to_cache(self, idx, trajectory_idx, replica_idx, bias):
        obj = Sample(
                trajectory=self.storage.trajectories[int(trajectory_idx)],
                replica=replica_idx,
                bias=bias
            )
        obj.idx[self.storage] = idx
        obj._origin = self.storage

        del obj.details
        del obj.ensemble
        del obj.mover
        del obj.parent

        self.cache[idx] = obj

        return obj


class SampleSetStore(ObjectStore):

    def __init__(self, storage):
        super(SampleSetStore, self).__init__(storage, SampleSet, json=False, load_partial=True)

        self.set_variable_partial_loading('movepath', self.update_movepath)

    update_movepath = func_update_object('movepath', 'sampleset', 'movepath', 'pathmovechanges')

    def save(self, sample_set, idx=None):
        # Check if all samples are saved
        map(self.storage.samples.save, sample_set)

        values = self.list_to_numpy(sample_set, 'sample')
        self.storage.variables['sampleset_sample_idx'][idx] = values

        self.storage.pathmovechanges.save(sample_set.movepath)
        self.save_object('sampleset_movepath', idx, sample_set.movepath)

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
        samples = self.list_from_numpy(values, 'samples')
        sample_set = SampleSet(samples, movepath=self.storage.pathmovechanges.load(movepath_idx))

        return sample_set

    def load_empty(self, idx):
        values = self.storage.variables['sampleset_sample_idx'][idx]

        # typecast to sample
        samples = self.list_from_numpy(values, 'samples')
        sample_set = SampleSet(samples, movepath=None)

        del sample_set.movepath

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

    def cache_all(self):
        """Load all samples as fast as possible into the cache

        """
        if not self._cached_all:
            idxs = range(len(self))
            values = self.storage.variables['sampleset_sample_idx'][:]

            # assume that these are cached!
            all_samples = self.storage.samples


            [ self.add_empty_to_cache(i,t,all_samples) for i,t in zip(
                idxs,
                values
                ) ]

            self._cached_all = True


    def add_empty_to_cache(self, idx, sample_idxs, all_samples):
        if idx not in self.cache:
            obj = SampleSet(
                    samples=[all_samples[sample_idx.tolist()] for sample_idx in sample_idxs],
                    movepath=None
                )
            obj.idx[self.storage] = idx
            obj._origin = self.storage

            del obj.movepath

            self.cache[idx] = obj