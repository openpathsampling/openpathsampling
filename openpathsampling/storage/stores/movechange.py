from openpathsampling.movechange import MoveChange
from openpathsampling.netcdfplus import StorableObject, ObjectStore

from uuid import UUID


class MoveChangeStore(ObjectStore):
    def __init__(self):
        super(MoveChangeStore, self).__init__(
            MoveChange,
            json=False
        )

        self._cached_all = False
        self.class_list = StorableObject.objects()

    def to_dict(self):
        return {}

    def _save(self, movechange, idx):
        self.vars['samples'][idx] = movechange.samples
        self.vars['input_samples'][idx] = movechange.input_samples
        self.vars['subchanges'][idx] = movechange.subchanges
        self.write('details', idx, movechange)
        self.vars['mover'][idx] = movechange.mover
        self.vars['cls'][idx] = movechange.__class__.__name__

    def _load(self, idx):
        cls_name = self.vars['cls'][idx]

        cls = self.class_list[cls_name]
        obj = cls.__new__(cls)
        MoveChange.__init__(obj, mover=self.vars['mover'][idx])

        obj.samples = self.vars['samples'][idx]
        obj.subchanges = self.vars['subchanges'][idx]
        obj.details = self.vars['details'][idx]
        try:
            obj.input_samples = self.vars['input_samples'][idx]
        except KeyError:  # BACKWARDS COMPATIBILITY; REMOVE IN 2.0
            obj.input_samples = None

        return obj

    def initialize(self, units=None):
        super(MoveChangeStore, self).initialize()

        # New short-hand definition
        self.create_variable('details', 'lazyobj.details')
        self.create_variable('mover', 'obj.pathmovers')
        self.create_variable('cls', 'str')

        self.create_variable('subchanges', 'obj.movechanges',
                             dimensions='...',
                             chunksizes=(65536,))

        self.create_variable('samples', 'obj.samples',
                             dimensions='...',
                             chunksizes=(65536,))

        self.create_variable('input_samples', 'obj.samples',
                             dimensions='...',
                             chunksizes=(10240,))

    def cache_all(self):
        """Load all samples as fast as possible into the cache

        """
        if not self._cached_all:
            poss = range(len(self))
            uuids = self.vars['uuid']

            cls_names = self.variables['cls'][:]
            samples_idxss = self.variables['samples'][:]
            subchanges_idxss = self.variables['subchanges'][:]
            mover_idxs = self.variables['mover'][:]
            details_idxs = self.variables['details'][:]
            try:
                input_samples_vars = self.variables['input_samples']
            except KeyError:
                # BACKWARD COMPATIBILITY: REMOVE IN 2.0
                input_samples_idxss = [[] for _ in samples_idxss]
            else:
                input_samples_idxss = input_samples_vars[:]

            [self._add_empty_to_cache(*v) for v in zip(
                poss,
                uuids,
                cls_names,
                samples_idxss,
                input_samples_idxss,
                mover_idxs,
                details_idxs)]

            [self._load_partial_subchanges(c, s) for c, s in zip(
                self,
                subchanges_idxss)]

            self._cached_all = True

    def _add_empty_to_cache(self, pos, uuid, cls_name, samples_idxs,
                            input_samples_idxs, mover_idx, details_idx):

        if pos not in self.cache:
            obj = self._load_partial_samples(cls_name, samples_idxs,
                                             input_samples_idxs, mover_idx,
                                             details_idx)

            obj.__uuid__ = uuid
            self.cache[pos] = obj

    def _load_partial_subchanges(self, obj, subchanges_idxs):
        if len(subchanges_idxs) > 0:
            subchanges_idxs = self.storage.to_uuid_chunks(subchanges_idxs)
            obj.subchanges = \
                [self.load(int(UUID(idx))) for idx in subchanges_idxs]

        return obj

    def _load_partial_samples(self, cls_name, samples_idxs,
                              input_samples_idxs, mover_idx, details_idx):
        cls = self.class_list[cls_name]
        obj = cls.__new__(cls)
        MoveChange.__init__(obj)

        if mover_idx[0] != '-':
            obj.mover = self.storage.pathmovers.load(int(UUID(mover_idx)))

        if len(samples_idxs) > 0:
            samples_idxs = self.storage.to_uuid_chunks(samples_idxs)
            obj.samples = [
                self.storage.samples.load(int(UUID(idx)))
                for idx in samples_idxs]
        else:
            obj.samples = []

        if len(input_samples_idxs) > 0:
            input_samples_idxs = \
                self.storage.to_uuid_chunks(input_samples_idxs)
            obj.input_samples = [
                self.storage.samples.load(int(UUID(idx))) if idx[0] != '-' else
                None for idx in input_samples_idxs]
        else:
            obj.input_samples = []

        if details_idx[0] != '-':
            obj.details = self.storage.details.proxy(int(UUID(details_idx)))

        return obj
