from openpathsampling.storage import ObjectStore
from openpathsampling.pathsimulator import MCStep


class MCStepStore(ObjectStore):
    def __init__(self):
        super(MCStepStore, self).__init__(
            MCStep,
            json=False
        )

        self._cached_all = False

    def _save(self, mcstep, idx):
        self.vars['change'][idx] = mcstep.change
        self.vars['active'][idx] = mcstep.active
        self.vars['previous'][idx] = mcstep.previous
        self.vars['simulation'][idx] = mcstep.simulation
        self.vars['mccycle'][idx] = mcstep.mccycle

    def _load(self, idx):
        return MCStep(
            mccycle=self.vars['mccycle'][idx],
            previous=self.vars['previous'][idx],
            active=self.vars['active'][idx],
            simulation=self.vars['simulation'][idx],
            change=self.vars['change'][idx]
        )

    def _init(self, units=None):
        super(MCStepStore, self)._init()

        # New short-hand definition
        self.init_variable('change', 'obj.pathmovechanges', chunksizes=(1,))
        self.init_variable('active', 'obj.samplesets', chunksizes=(1,))
        self.init_variable('previous', 'obj.samplesets', chunksizes=(1,))
        self.init_variable('simulation', 'obj.pathsimulators', chunksizes=(1,))
        self.init_variable('mccycle', 'int', chunksizes=(1,))

    def all(self):
        self.cache_all()
        return self

    def cache_all(self):
        """Load all samples as fast as possible into the cache

        """
        if not self._cached_all:
            idxs = range(len(self))

            storage = self.storage

            steps = storage.variables[self.prefix + '_mccycle'][:]
            previous_idxs = storage.variables[self.prefix + '_previous'][:]
            active_idxs = storage.variables[self.prefix + '_active'][:]
            simulation_idxs = storage.variables[self.prefix + '_simulation'][:]
            change_idxs = storage.variables[self.prefix + '_change'][:]

            [self.add_to_cache(*v) for v in zip(
                idxs,
                steps,
                previous_idxs,
                active_idxs,
                simulation_idxs,
                change_idxs)]

            self._cached_all = True

    def add_to_cache(self, idx, step, previous_idx,
                     active_idx, simulation_idx, change_idx):
        if idx not in self.cache:
            storage = self.storage
            obj = MCStep(
                mccycle=int(step),
                previous=storage.samplesets[int(previous_idx)],
                active=storage.samplesets[int(active_idx)],
                simulation=storage.pathsimulators[int(simulation_idx)],
                change=storage.pathmovechanges[int(change_idx)]
            )

            self.index[obj] = idx
            self.cache[idx] = obj
