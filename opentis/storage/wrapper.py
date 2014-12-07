MULTIFILE = False

def savecache(func):
    def inner(self, obj, idx = None, *args, **kwargs):

        if MULTIFILE:
            # Check, if in multifile this is the right storage, otherwise skip to the next or previous one.

            if self._min_idx > idx:
                # use the previous storage
                store = self.storage._previous._storages[self.__class__.__name__]
                store.save(obj, idx, *args, **kwargs)
            elif self._max_idx is not None and self._max_idx < idx:
                store = self.storage._next._storages[self.__class__.__name__]
                store.save(obj, idx, *args, **kwargs)


        idx = self.index(obj, idx)
        if idx is not None:
            func(self, obj, idx, *args, **kwargs)

        # store the ID in the cache
        self.cache[idx] = obj
        if self.named and hasattr(obj, 'name') and obj.name != '':
            self.cache[obj.name] = obj

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
                self.all_names[obj.identifier] = find_idx
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
        # TODO: Maybe this functionality should be in a separate function

        if type(idx) is not str and idx < 0:
            return None


        if MULTIFILE:
            # Check, if in multifile this is the right storage, otherwise skip to the next or previous one.

            if self._min_idx > idx:
                # use the previous storage
                store = self.storage._previous._storages[self.__class__.__name__]
                return store.load(idx, *args, **kwargs)
            elif self._max_idx is not None and self._max_idx < idx:
                store = self.storage._next._storages[self.__class__.__name__]
                store.load(idx, *args, **kwargs)

        n_idx = idx

        if idx in self.cache:
            cc = self.cache[idx]
            if type(cc) is int:
                # this happens when we want to load by name (str)
                # and we need to actually load it
                n_idx = cc
            else:

                return self.cache[idx]
        elif type(idx) is str:
            if self.named:
                if not self._names_loaded:
                    self.update_name_cache()
                if idx in self.cache:
                    n_idx = self.cache[idx]
                else:
                    raise ValueError('str "' + idx + '" not found in storage')
            else:
                raise ValueError('str "' + idx + '" as indices are only allowed in named storage')

            n_idx = self.idx_from_name(idx)

        if MULTIFILE:
            obj = func(self, n_idx - self._min_idx, *args, **kwargs)
        else:
            obj = func(self, n_idx, *args, **kwargs)

        obj.idx[self.storage] = n_idx

        self.cache[obj.idx[self.storage]] = obj


        if self.named and hasattr(obj, 'name') and obj.name != '':
            self.cache[obj.name] = obj

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