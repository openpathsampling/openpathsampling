from openpathsampling.storage import ObjectStore
from openpathsampling.pathmovechange import PathMoveChange
from openpathsampling.storage.object_storage import func_update_object
from openpathsampling.todict import class_list

class PathMoveChangeStore(ObjectStore):
    def __init__(self, storage):
        super(PathMoveChangeStore, self).__init__(
            storage,
            PathMoveChange,
            json=False,
            load_partial=True
        )

        self.set_variable_partial_loading('details', self.update_details)
#        self.set_variable_partial_loading('mover', self.update_mover)

        self._cached_all = False

    def load_empty(self, idx):

        obj = self._load_partial(idx)
        del obj.details

        return obj

    update_details = func_update_object('details', 'change', 'details', '_details')
    update_mover = func_update_object('mover', 'change', 'pathmover', 'pathmovers')


    def save(self, pathmovechange, idx=None):
        if idx is not None:
            if len(pathmovechange.samples) > 0:
                map(self.storage.samples.save, pathmovechange.samples)

            values = self.list_to_numpy(pathmovechange.samples, 'samples')
            self.storage.variables['change_samples_idxs'][idx] = values

            if len(pathmovechange.subchanges) > 0:
                map(self.storage.pathmovechanges.save, pathmovechange.subchanges)

            values = self.list_to_numpy(pathmovechange.subchanges, 'pathmovechanges')
            self.storage.variables['change_subchanges_idxs'][idx] = values

            self.save_object('change_details', idx, pathmovechange.details)
            self.save_object('change_pathmover', idx, pathmovechange.mover)

            self.save_variable('change_cls', idx, pathmovechange.__class__.__name__)

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

        obj = self._load_partial(idx)

        details_idx = int(self.storage.variables['change_details_idx'][idx])
        obj.details = self.storage._details[details_idx]

        return obj

    def _init(self, units=None):
        super(PathMoveChangeStore, self)._init()

        # New short-hand definition
        self.init_variable('change_details_idx', 'index', chunksizes=(1, ))
        self.init_variable('change_pathmover_idx', 'index', chunksizes=(1, ))
        self.init_variable('change_cls', 'str', chunksizes=(1, ))

        self.init_variable('change_subchanges_idxs', 'index',
            variable_length = True,
            chunksizes=(10240, )
        )

        self.init_variable('change_samples_idxs', 'index',
            variable_length = True,
            chunksizes=(10240, )
        )

    def all(self):
        self.cache_all()
        return self

    def cache_all(self):
        """Load all samples as fast as possible into the cache

        """
        if not self._cached_all:
            idxs = range(len(self))

            cls_names = self.storage.variables['change_cls'][:]
            samples_idxss = self.storage.variables['change_samples_idxs'][:]
            subchanges_idxss = self.storage.variables['change_subchanges_idxs'][:]
            mover_idxs = self.storage.variables['change_pathmover_idx'][:]

            [ self.add_empty_to_cache(i,c,g,m) for i,c,g,m in zip(
                idxs,
                cls_names,
                samples_idxss,
                mover_idxs) ]

            [ self._load_partial_subchanges(c,s) for c,s in zip(
                self,
                subchanges_idxss) ]

            self._cached_all = True


    def add_empty_to_cache(self, idx, cls_name, samples_idxs, mover_idx):

        if idx not in self.cache:
            obj = self._load_partial_samples(cls_name, samples_idxs, mover_idx)
            obj.idx[self.storage] = idx
            obj._origin = self.storage

            self.cache[idx] = obj
            del obj.details
#            del obj.mover



    def _load_partial(self, idx):
        samples_idxs = self.storage.variables['change_samples_idxs'][idx]
        subchanges_idxs = self.storage.variables['change_subchanges_idxs'][idx]
        mover_idx = self.storage.variables['change_pathmover_idx'][idx]

        cls_name = self.storage.variables['change_cls'][idx]

        obj = self._load_partial_samples(cls_name, samples_idxs, mover_idx)
        return self._load_partial_subchanges(obj, subchanges_idxs)

    def _load_partial_subchanges(self, obj, subchanges_idxs):
        if len(subchanges_idxs) > 0:
            obj.subchanges = [ self.load(int(idx)) for idx in subchanges_idxs ]

        return obj

    def _load_partial_samples(self, cls_name, samples_idxs, mover_idx):
        cls = class_list[cls_name]
        obj = cls.__new__(cls)
        PathMoveChange.__init__(obj, mover=self.storage.pathmovers[int(mover_idx)])

        if len(samples_idxs) > 0:
            obj.samples = [ self.storage.samples[int(idx)] for idx in samples_idxs ]

        return obj
