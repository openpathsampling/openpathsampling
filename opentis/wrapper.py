def storable(super_class):
    super_class.cls = super_class.__name__.lower()
    super_class.default_storage = None

    def _init(self, *args, **kwargs):
        super_class._init(self, *args, **kwargs)

        if 'idx' in kwargs:
            self.idx = kwargs['idx']
        else:
            self.idx = dict()

        if 'storage' in kwargs:
            self.default_storage = kwargs['storage']

    def _save(self, storage=None):
        if storage is None:
            storage = self.default_storage
        storage.save(self)

    super_class._init = super_class.__init__
    super_class.__init__ = _init

    super_class.save = _save

    return super_class