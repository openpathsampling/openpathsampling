from openpathsampling.storage import ObjectStore

class NameStore(ObjectStore):
    def __init__(self):
        super(NameStore, self).__init__(
            type(None),
            json=False
        )

        self.names = dict()
        self.stores = dict()

    def save(self, obj, idx=None, name=None):
        if name is None:
            name = obj.name

        obj.save(self.storage)
        self.vars['name'][idx] = name
        obj.fix_name()
        store = self.storage.find_store(obj)
        obj_idx = self.storage.index.get(obj, None)
        self.vars['store'][idx] = store
        self.vars['idx'][idx] = obj_idx

        self.names[idx] = name
        self.stores[idx] = store

    def load(self, idx):
        store = self.vars['store'][idx],
        idx = self.vars['idx'][idx],

        return store.load(idx)

    def find(self, needle, store=None):
        if store is None:
            return [self.load(idx) for idx, name in self.names.iteritems() if name == needle]
        else:
            return [self.load(idx) for idx in self.names.keys()
                    if self.names[idx] == needle and self.stores[idx] == store]

    def _restore(self):
        """
        This will cache all names
        """
        self.names = dict(zip(*[range(len(self)) , self.vars['name'][:]]))
        self.stores = dict(zip(*[range(len(self)) , self.vars['store'][:]]))

    def _init(self, units=None):
        super(NameStore, self)._init()

        # New short-hand definition
        self.init_variable('name', 'str', chunksizes=(10240, ))
        self.init_variable('store', 'store', chunksizes=(10240, ))
        self.init_variable('idx', 'index', chunksizes=(100, ))

    def __getitem__(self, item):
        if type(item) is str:
            return self.find(item)