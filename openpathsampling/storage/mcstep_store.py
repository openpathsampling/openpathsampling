from openpathsampling.storage import ObjectStore
from openpathsampling.pathsimulator import MCStep

class MCStepStore(ObjectStore):
    def __init__(self, storage):
        super(MCStepStore, self).__init__(
            storage,
            MCStep,
            json=False,
            load_partial=False
        )

        self._cached_all = False

    def save(self, mcstep, idx=None):
        if idx is not None:
            self.save_object('mcstep_change', idx, mcstep.change)
            self.save_object('mcstep_active', idx, mcstep.active)
            self.save_object('mcstep_previous', idx, mcstep.previous)
            self.save_object('mcstep_simulation', idx, mcstep.simulation)
            self.save_object('mcstep_step', idx, mcstep.step_number)

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

        storage = self.storage

        previous_idx = self.storage.variables['mcstep_previous_idx'][idx]
        active_idx = self.storage.variables['mcstep_active_idx'][idx]
        simulation_idx = self.storage.variables['mcstep_simulation_idx'][idx]
        change_idx = self.storage.variables['mcstep_simulation_idx'][idx]

        step = self.storage.variables['mcstep_step'][idx]

        return MCStep(
            step_number=step,
            previous=storage.samplesets[previous_idx],
            active=storage.samplesets[active_idx],
            simulation=storage.samplesets[simulation_idx],
            change=storage.samplesets[change_idx]
        )

    def _init(self, units=None):
        super(MCStepStore, self)._init()

        # New short-hand definition
        self.init_variable('mcstep_change_idx', 'index', chunksizes=(1, ))
        self.init_variable('mcstep_active_idx', 'index', chunksizes=(1, ))
        self.init_variable('mcstep_previous_idx', 'index', chunksizes=(1, ))
        self.init_variable('mcstep_simulation_idx', 'index', chunksizes=(1, ))
        self.init_variable('mcstep_step', 'int', chunksizes=(1, ))

    def all(self):
        self.cache_all()
        return self

    def cache_all(self):
        """Load all samples as fast as possible into the cache

        """
        if not self._cached_all:
            idxs = range(len(self))

            storage = self.storage

            steps = storage.variables['change_cls'][:]
            previous_idxs = storage.variables['change_samples_idxs'][:]
            active_idxs = storage.variables['change_subchanges_idxs'][:]
            simulation_idxs = storage.variables['change_pathmover_idx'][:]
            change_idxs = storage.variables['']

            [ self.add_to_cache(i, n, p, a, s, c) for i, n, p, a, s, c in zip(
                idxs,
                steps,
                previous_idxs,
                active_idxs,
                simulation_idxs,
                change_idxs) ]

            self._cached_all = True


    def add_to_cache(self, idx, step, previous_idx,
                     active_idx, simulation_idx, change_idx):
        if idx not in self.cache:

            storage = self.storage
            obj = MCStep(
                step_number=step,
                previous=storage.samplesets[previous_idx],
                active=storage.samplesets[active_idx],
                simulation=storage.samplesets[simulation_idx],
                change=storage.samplesets[change_idx]
            )
            obj.idx[self.storage] = idx

            self.cache[idx] = obj