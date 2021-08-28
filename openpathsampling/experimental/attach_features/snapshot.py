# TODO: get these from some universal source
def set_uuid(obj, uuid):
    obj.__uuid__ = uuid

def get_uuid(obj):
    return obj.__uuid__


class Snapshot(StorableObject):
    def __init__(self, **kwargs):
        super().__init__()
        self._reversed_uuid = self.reverse_uuid()  # TODO: change SimStore
        self._reversed_snapshot = None
        for arg, val in kwargs.items():
            setattr(self, arg, val)

    @classmethod
    def _is_compatible(cls, other):
        """check that all required features of ``self`` are in ``other``"""
        required = [p.name for p in cls.__features__.parameters
                    if p.default is p.empty]
        for name in required:
            try:
                _ = getattr(other, name)
            except AttributeError:
                return False

        return True

    def copy_with_replacement(self, **kwargs):
        dct = self.to_dict()

        # copy things that need to be copied
        copy_params = [p for p in self.__features__.parameters
                       if p not in kwargs]
        for param in copy_params:
            dct[param.name] = param.copy(dct[param.name])

        # update with new versions
        dct.update(kwargs)
        return self.from_dict(dct)

    def to_dict(self):
        return {p.name: getattr(self, p.name)
                for p in self.__features__.parameters}

    def _create_reversed(self):
        reversed_dct = {p.name: p.time_reverse(getattr(self, p.name))
                        for p in self.__features__.parameters
                        if 'time_reverse' in p.operations}
        rev = self.copy_with_replacement(**reversed_dct)
        set_uuid(rev, self._reversed_uuid)
        rev._reversed_uuid = get_uuid(self)
        rev._reversed_snapshot = self
        self._reversed_snapshot = rev
        return rev

    @property
    def reverse(self):
        if self._reversed_snapshot is None:
            self._reversed_snapshot = self._create_reversed()
        return self._reversed_snapshot

    @classmethod
    def create_empty(cls):
        new = cls.__new__(cls)
        super(cls, new).__init__()


