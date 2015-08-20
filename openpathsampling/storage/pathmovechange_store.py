from openpathsampling.storage import ObjectStore
from openpathsampling.pathmovechange import PathMoveChange
from openpathsampling.todict import OPSObject

class PathMoveChangeStore(ObjectStore):
    def __init__(self):
        super(PathMoveChangeStore, self).__init__(
            PathMoveChange,
            json=False,
            load_partial=True
        )

        self.set_variable_partial_loading('details')
#        self.set_variable_partial_loading('mover')

        self._cached_all = False
        self.class_list = OPSObject.objects()

    def load_empty(self, idx):

        obj = self._load_partial(idx)
        del obj.details

        return obj

    def save(self, pathmovechange, idx=None):
        if idx is not None:
            self.vars['samples'][idx] = pathmovechange.samples
            self.vars['subchanges'][idx] = pathmovechange.subchanges
            self.vars['details'][idx] = pathmovechange.details
            self.vars['mover'][idx] = pathmovechange.mover
            self.vars['cls'][idx] = pathmovechange.__class__.__name__

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
        obj.details = self.vars['details'][idx]

        return obj

    def _init(self, units=None):
        super(PathMoveChangeStore, self)._init()

        # New short-hand definition
        self.init_variable('details', 'obj.details', chunksizes=(1, ))
        self.init_variable('mover', 'obj.pathmovers', chunksizes=(1, ))
        self.init_variable('cls', 'str', chunksizes=(1, ))

        self.init_variable('subchanges', 'obj.pathmovechanges',
            variable_length = True,
            chunksizes=(10240, )
        )

        self.init_variable('samples', 'obj.samples',
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

            cls_names = self.storage.variables[self.prefix + '_cls'][:]
            samples_idxss = self.storage.variables[self.prefix + '_samples'][:]
            subchanges_idxss = self.storage.variables[self.prefix + '_subchanges'][:]
            mover_idxs = self.storage.variables[self.prefix + '_mover'][:]

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
            obj.idx[self] = idx
            obj._origin = self

            self.cache[idx] = obj
            del obj.details
#            del obj.mover

    def _load_partial(self, idx):
        cls_name = self.vars['cls'][idx]

        cls = self.class_list[cls_name]
        obj = cls.__new__(cls)
        PathMoveChange.__init__(obj, mover=self.vars['mover'][idx])

        obj.samples = self.vars['samples'][idx]
        obj.subchanges = self.vars['subchanges'][idx]

        return obj

    def _load_partial_subchanges(self, obj, subchanges_idxs):
        if len(subchanges_idxs) > 0:
            obj.subchanges = [ self.load(int(idx)) for idx in subchanges_idxs ]

        return obj

    def _load_partial_samples(self, cls_name, samples_idxs, mover_idx):
        cls = self.class_list[cls_name]
        obj = cls.__new__(cls)
        PathMoveChange.__init__(obj, mover=self.storage.pathmovers[int(mover_idx)])

        if len(samples_idxs) > 0:
            obj.samples = [ self.storage.samples[int(idx)] for idx in samples_idxs ]

        return obj
