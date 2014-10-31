def storable(super_class):
    class NewClass(super_class):
        cls = super_class.__name__.lower()
        default_storage = None

        def __init__(self, *args, **kwargs):
            super(NewClass, self).__init__(*args, **kwargs)
            if 'idx' in kwargs:
                self.idx = kwargs['idx']
            else:
                self.idx = dict()

        def save(self, storage=None):
            if storage is None:
                storage = self.default_storage
            storage.save(self)

    NewClass.__name__ = super_class.__name__
    return NewClass