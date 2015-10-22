from object_storage import ObjectStore
from openpathsampling.sample import SampleSet, Sample
from objproxy import LoaderProxy


class SampleStore(ObjectStore):
    def __init__(self):
        super(SampleStore, self).__init__(Sample, json=False)
        self._cached_all = False

    def _save(self, sample, idx):
        self.vars['trajectory'][idx] = sample.trajectory
        self.vars['ensemble'][idx] = sample.ensemble
        self.vars['replica'][idx] = sample.replica
        self.write('parent', idx, sample)
        self.write('details', idx, sample)
        self.vars['bias'][idx] = sample.bias
        self.vars['mover'][idx] = sample.mover

    def _load(self, idx):
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
        return [sample for sample in self.iterator() if sample.ensemble == ensemble]

    def _init(self):
        super(SampleStore, self)._init()

        # New short-hand definition
        self.init_variable('trajectory', 'obj.trajectories', chunksizes=(1,))
        self.init_variable('ensemble', 'obj.ensembles', chunksizes=(1,))
        self.init_variable('replica', 'int', chunksizes=(1,))
        self.init_variable('parent', 'lazyobj.samples', chunksizes=(1,))
        self.init_variable('details', 'lazyobj.details', chunksizes=(1,))
        self.init_variable('bias', 'float', chunksizes=(1,))
        self.init_variable('mover', 'obj.pathmovers', chunksizes=(1,))

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

            [self._add_empty_to_cache(*v) for v in zip(
                idxs,
                trajectory_idxs,
                replica_idxs,
                biass,
                ensemble_idxs,
                parent_idxs,
                details_idxs,
                mover_idxs)]

            self._cached_all = True

    def _add_empty_to_cache(self, idx, trajectory_idx, replica_idx, bias,
                           ensemble_idx, parent_idx, details_idx, mover_idx):
        storage = self.storage
        obj = Sample(
            trajectory=storage.trajectories[int(trajectory_idx)],
            replica=int(replica_idx),
            bias=float(bias),
            ensemble=storage.ensembles[int(ensemble_idx)],
            mover=storage.pathmovers[int(mover_idx)],
            parent=LoaderProxy(storage.samples, int(parent_idx)),
            details=LoaderProxy(storage.details, int(details_idx))
        )

        self.index[obj] = idx
        self.cache[idx] = obj

        return obj


class SampleSetStore(ObjectStore):
    def __init__(self):
        super(SampleSetStore, self).__init__(SampleSet, json=False)

    def _save(self, sample_set, idx):
        map(self.storage.samples.save, sample_set)

        self.vars['samples'][idx] = sample_set
        self.write('movepath', idx, sample_set)

    def _load(self, idx):
        sample_set = SampleSet(
            self.vars['samples'][idx],
            movepath=LoaderProxy(self.storage.pathmovechanges, int(self.variables['movepath'][idx]))
        )

        return sample_set

    def sample_indices(self, idx):
        """
        Load sample indices for sample_set with ID 'idx' from the storage

        Parameters
        ----------
        idx : int
            ID of the sample_set

        Returns
        -------
        list of int
            list of sample indices
        """

        return self.variables['samples'][idx].tolist()

    def _init(self, units=None):
        """
        Initialize the associated storage to allow for sampleset storage

        """
        super(SampleSetStore, self)._init()

        self.init_variable('samples', 'obj.samples',
                           dimensions='...',
                           description="sampleset[sampleset][frame] is the sample index " +
                                       "(0..nspanshots-1) of frame 'frame' of sampleset 'sampleset'.",
                           chunksizes=(1024,)
                           )

        self.init_variable('movepath', 'obj.pathmovechanges', chunksizes=(1,))

    def cache_all(self):
        """Load all samples as fast as possible into the cache

        """
        if not self._cached_all:
            idxs = range(len(self))
            samples_idxs = self.variables['samples'][:]
            pmc_idxs = self.variables['movepath'][:]

            [self._add_empty_to_cache(*v) for v in zip(
                idxs,
                samples_idxs,
                pmc_idxs
            )]

            self._cached_all = True

    def _add_empty_to_cache(self, idx, sample_idxs, pmc_idx):
        if idx not in self.cache:
            obj = SampleSet(
                samples=[self.storage.samples[sample_idx.tolist()] for sample_idx in sample_idxs],
                movepath=LoaderProxy(self.storage.pathmovechanges, int(pmc_idx))
            )

            self.index[obj] = idx
            self.cache[idx] = obj
