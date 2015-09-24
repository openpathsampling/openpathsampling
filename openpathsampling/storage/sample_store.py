from object_storage import ObjectStore
from openpathsampling.sample import SampleSet, Sample
from objproxy import LoaderProxy

class SampleStore(ObjectStore):
    def __init__(self):
        super(SampleStore, self).__init__(Sample, json=False, load_partial=False)

#        self.set_variable_partial_loading('details')
#        self.set_variable_partial_loading('parent')
#        self.set_variable_partial_loading('mover')
#        self.set_variable_partial_loading('ensemble')

        self._cached_all = False


    def load_empty(self, idx):
        obj = Sample(
            trajectory=self.vars['trajectory'][idx],
            replica=self.vars['replica'][idx],
            bias=self.vars['bias'][idx],
        )

#        del obj.details
#        del obj.ensemble
#        del obj.mover
#        del obj.parent

        return obj

    def save(self, sample, idx=None):
        self.vars['trajectory'][idx] = sample.trajectory
        self.vars['ensemble'][idx] = sample.ensemble
        self.vars['replica'][idx] = sample.replica
        self.vars['parent'][idx] = sample.parent
        self.vars['details'][idx] = sample.details
        self.vars['bias'][idx] = sample.bias
        self.vars['mover'][idx] = sample.mover

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

        obj = Sample(
            trajectory=self.vars['trajectory'][idx],
            replica=self.vars['replica'][idx],
            ensemble=self.vars['ensemble'][idx],
            parent=self.vars['parent'][idx],
            details=self.vars['details'][idx],
            bias=self.vars['bias'][idx],
            mover=self.vars['mover'][idx]
        )

        return obj

    def by_ensemble(self, ensemble):
        return [ sample for sample in self.iterator() if sample.ensemble == ensemble ]

    def _init(self):
        super(SampleStore, self)._init()

        # New short-hand definition
        self.init_variable('trajectory', 'obj.trajectories', chunksizes=(1, ))
        self.init_variable('ensemble', 'obj.ensembles', chunksizes=(1, ))
        self.init_variable('replica', 'int', chunksizes=(1, ))
        self.init_variable('parent', 'lazyobj.samples', chunksizes=(1, ))
        self.init_variable('details', 'lazyobj.details', chunksizes=(1, ))
        self.init_variable('bias', 'float', chunksizes=(1, ))
        self.init_variable('mover', 'obj.pathmovers', chunksizes=(1, ))

    def all(self):
        self.cache_all()
        return self

    def cache_all(self):
        """Load all samples as fast as possible into the cache

        """
        if not self._cached_all:
            idxs = range(len(self))
            trajectory_idxs = self.variables['trajectory'][:]
            replica_idxs = self.variables['replica'][:]
            biass = self.variables['bias'][:]
            ensemble_idxs = self.variables['ensemble'][:]
            parent_idxs = self.variables['parent'][:]
            mover_idxs = self.variables['mover'][:]
            details_idxs = self.variables['details'][:]

            [ self.add_empty_to_cache(*v) for v in zip(
                idxs,
                trajectory_idxs,
                replica_idxs,
                biass,
                ensemble_idxs,
                parent_idxs,
                details_idxs,
                mover_idxs) ]

            self._cached_all = True

    def add_empty_to_cache(self, idx, trajectory_idx, replica_idx, bias,
                           ensemble_idx, parent_idx, details_idx, mover_idx):
        obj = Sample(
                trajectory=self.storage.trajectories[int(trajectory_idx)],
                replica=int(replica_idx),
                bias=float(bias),
                ensemble=self.storage.ensembles[int(ensemble_idx)],
                mover=self.storage.pathmovers[int(mover_idx)],
                parent=LoaderProxy({self.storage.samples: int(parent_idx)}),
                details=LoaderProxy({self.storage.details: int(details_idx)})
            )

        self.index[obj] = idx
        self.cache[idx] = obj

        return obj


class SampleSetStore(ObjectStore):

    def __init__(self):
        super(SampleSetStore, self).__init__(SampleSet, json=False, load_partial=True)

    def save(self, sample_set, idx=None):
        # Check if all samples are saved
        map(self.storage.samples.save, sample_set)

        self.vars['samples'][idx] = sample_set
        self.vars['movepath'][idx] = sample_set.movepath

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
        values = self.storage.variables[self.prefix + '_samples'][idx]

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

        sample_set = SampleSet(
            self.vars['samples'][idx],
            movepath=LoaderProxy({self.storage.pathmovechanges: int(self.variables['movepath'])})
        )

        return sample_set

    def _init(self, units=None):
        """
        Initialize the associated storage to allow for sampleset storage

        """
        super(SampleSetStore, self)._init()

        self.init_variable('samples', 'obj.samples',
            dimensions='...',
            description="sampleset[sampleset][frame] is the sample index (0..nspanshots-1) of frame 'frame' of sampleset 'sampleset'.",
            chunksizes=(1024, )
        )

        self.init_variable('movepath', 'obj.pathmovechanges', chunksizes=(1, ))

    def cache_all(self):
        """Load all samples as fast as possible into the cache

        """
        if not self._cached_all:
            idxs = range(len(self))
            samples_idxs = self.variables['samples'][:]
            pmc_idxs = self.variables['movepath'][:]

            [ self.add_empty_to_cache(*v) for v in zip(
                idxs,
                samples_idxs,
                pmc_idxs
                ) ]

            self._cached_all = True

    def add_empty_to_cache(self, idx, sample_idxs, pmc_idx):
        if idx not in self.cache:
            obj = SampleSet(
                    samples=[self.storage.samples[sample_idx.tolist()] for sample_idx in sample_idxs],
                    movepath=LoaderProxy({self.storage.pathmovechanges : int(pmc_idx)})
                )

            self.index[obj] = idx
            self.cache[idx] = obj