def storable(super_class):
    class NewClass(super_class):
        cls = super_class.__name__.lower()
        def __init__(self, *args, **kwargs):
            super(NewClass, self).__init__(*args, **kwargs)
            if 'idx' in kwargs:
                self.idx = kwargs['idx']
            else:
                self.idx = dict()

    NewClass.__name__ = super_class.__name__
    return NewClass


def savecache(func):
    def inner(self, obj, idx = None, *args, **kwargs):
        idx = self.index(obj, idx)
        if idx is not None:
            func(self, obj, idx, *args, **kwargs)
    return inner


def identifiable(func):
    def inner(self, obj, idx=None, *args, **kwargs):
        if idx is None and hasattr(obj, 'identifier'):
            if not hasattr(obj,'json'):
                setattr(obj,'json',self.get_json(obj))

            find_idx = self.find_by_identifier(obj.identifier)
            if find_idx is not None:
                # found and does not need to be saved, but we will let this ensemble point to the storage
                # in case we want to save and need the idx
                obj.idx[self.storage] = find_idx
                self.cache[find_idx] = obj
            else:
                func(self, obj, idx, *args, **kwargs)
                # Finally register with the new idx in the identifier cache dict.
                new_idx = obj.idx[self.storage]
                self.all_names[obj.identifier] = new_idx
        else:
            func(self, obj, idx, *args, **kwargs)

    return inner


def loadcache(func):
    def inner(self, idx, *args, **kwargs):
        if idx in self.cache:
            return self.cache[idx]

        obj = func(self, idx, *args, **kwargs)
        self.cache[idx] = obj
        return obj
    return inner