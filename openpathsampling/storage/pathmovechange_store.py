from openpathsampling.storage import ObjectStore
from openpathsampling.pathmovechange import PathMoveChange, SamplePathMoveChange
from openpathsampling.storage.object_storage import func_update_object
from openpathsampling.todict import class_list
import openpathsampling.pathmovechange as pmc

class PathMoveChangeStore(ObjectStore):
    def __init__(self, storage):
        super(PathMoveChangeStore, self).__init__(
            storage,
            PathMoveChange,
            json=False,
            load_partial=True
        )

        self.set_variable_partial_loading('details', self.update_details)
        self.set_variable_partial_loading('mover', self.update_mover)

        self._cached_all = False

    def load_empty(self, idx):

        obj = self._load_partial(idx)

        del obj.details
        del obj.mover

        return obj

    update_details = func_update_object('details', 'change', 'details', '_details')
    update_mover = func_update_object('mover', 'change', 'pathmover', 'pathmovers')


    def save(self, pathmovechange, idx=None):
        if idx is not None:
            if len(pathmovechange.generated) > 0:
                map(self.storage.samples.save, pathmovechange.generated)

            values = self.list_to_numpy(pathmovechange.generated, 'samples')
            self.storage.variables['change_generated_idxs'][idx] = values

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
        pathmover_idx = int(self.storage.variables['change_pathmover_idx'][idx])

        obj.mover = self.storage.pathmovers[pathmover_idx],
        obj.details = self.storage._details[details_idx]

        return obj

    def _init(self, units=None):
        super(PathMoveChangeStore, self)._init(units)

        # New short-hand definition
        self.init_variable('change_details_idx', 'index', chunksizes=(1, ))
        self.init_variable('change_pathmover_idx', 'index', chunksizes=(1, ))
        self.init_variable('change_cls', 'str', chunksizes=(1, ))

        self.init_variable('change_subchanges_idxs', 'index',
            variable_length = True,
            chunksizes=(10240, )
        )

        self.init_variable('change_generated_idxs', 'index',
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

            cls_names = self.load_variable('change_cls')[:]
            generated_idxss = int(self.storage.variables['change_generated_idxs'][:])
            subchanges_idxss = int(self.storage.variables['change_subchanges_idxs'][:])

            [ self.add_empty_to_cache(i,c,g,s) for i,c,g,s in zip(
                idxs,
                cls_names,
                generated_idxss,
                subchanges_idxss) ]

            self._cached_all = True


    def add_empty_to_cache(self, idx, cls_name, generated_idxs, subchanges_idxs):

        obj = self._load_partial(idx)

        del obj.details
        del obj.mover

        obj.idx[self.storage] = idx
        obj._origin = self.storage

        self.cache[idx] = obj

        return obj

    def _load_partial(self, idx):
        generated_idxs = self.storage.variables['change_generated_idxs'][idx]
        subchanges_idxs = self.storage.variables['change_subchanges_idxs'][idx]

        cls_name = self.storage.variables['change_cls'][idx]
        cls = class_list[cls_name]

#        obj = cls.__new__(cls)
#
#        if len(generated_idxs) > 0:
#            obj.generated = [ self.storage.samples[idx] for idx in generated_idxs ]
#
#        if len(subchanges_idxs) > 0:
#            obj.subchanges = [ self[idx] for idx in subchanges_idxs ]
#
#        print cls_name, subchanges_idxs

        if cls is pmc.SamplePathMoveChange:
            obj = cls(
                generated=[ self.storage.samples[int(idx)] for idx in generated_idxs ]
            )
        elif cls is pmc.EmptyPathMoveChange:
            obj = cls()
        elif issubclass(cls, pmc.SequentialPathMoveChange):
            obj = cls(
                subchanges=[ self.load(idx) for idx in subchanges_idxs ]
            )
        else:
            obj = cls(
                subchange=self.load(subchanges_idxs[0])
            )

        return obj

