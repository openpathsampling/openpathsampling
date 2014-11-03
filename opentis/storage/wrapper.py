def savecache(func):
    def inner(self, obj, idx = None, *args, **kwargs):
        idx = self.index(obj, idx)
        if idx is not None:
            func(self, obj, idx, *args, **kwargs)
    return inner


def saveidentifiable(func):
    def inner(self, obj, idx=None, *args, **kwargs):
        if idx is None and hasattr(obj, 'identifier'):
            if not hasattr(obj,'json'):
                setattr(obj,'json',self.object_to_json(obj))

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

def loadidentifiable(func):
    def inner(self, idx=None, *args, **kwargs):
        if idx is not None and type(idx) is str:
            find_idx = self.find_by_identifier(idx)
            if find_idx is not None:
                # names id is found so load with normal id
                return func(self, find_idx, *args, **kwargs)
            else:
                # named id does not exist
                return None
        else:
            return func(self, idx, *args, **kwargs)

    return inner
