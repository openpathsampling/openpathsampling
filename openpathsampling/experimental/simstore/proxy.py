from .uuids import get_uuid, set_uuid

import logging
logger = logging.getLogger(__name__)

def unproxy(uuid_mapping):
    """This is a convenience for unproxying a lot of objects at once.
    """
    for uuid, obj in uuid_mapping.items():
        if isinstance(obj, GenericLazyLoader):
            uuid_mapping[uuid] = obj.load()

def make_lazy_class(cls_):
    # this is to mix-in inheritence
    class LazyLoader(GenericLazyLoader, cls_):
        pass
    return LazyLoader

class GenericLazyLoader(object):
    def __init__(self, uuid, class_, storage):
        set_uuid(self, uuid)
        self.storage = storage
        self.class_ = class_
        self._loaded_object = None

    def load(self):
        if self._loaded_object is None:
            self._loaded_object = \
                    self.storage.load([get_uuid(self)], force=True)[0]
        if self._loaded_object is None:
            raise RuntimeError("UUID not found in storage: "
                               + get_uuid(self))
        return self._loaded_object

    def __getattr__(self, attr):
        # apparently IPython pretty-printing looks for a bunch of
        # attributes; this means we auto-load if we try to autoprint the
        # repr in IPython (TODO)
        return getattr(self.load(), attr)

    def __getitem__(self, item):
        return self.load()[item]

    def __iter__(self):
        return self.load().__iter__()

    def __len__(self):
        return len(self.load())

    def __str__(self):
        if self._loaded_object:
            return str(self._loaded_object)
        else:
            return repr(self)

    def __repr__(self):
        if self._loaded_object:
            return repr(self._loaded_object)
        else:
            return ("<LazyLoader for " + str(self.class_.__name__)
                    + " UUID " + str(self.__uuid__) + ">")


class ProxyObjectFactory(object):
    def __init__(self, storage, serialization_schema):
        self.storage = storage
        self.serialization_schema = serialization_schema
        self.lazy_classes = {}

    def make_lazy(self, cls, uuid):
        if cls not in self.lazy_classes:
            self.lazy_classes[cls] = make_lazy_class(cls)
        return self.lazy_classes[cls](uuid=uuid,
                                      class_=cls,
                                      storage=self.storage)

    def make_all_lazies(self, lazies):
        # lazies is dict of {table_name: list_of_lazy_uuid_rows}
        all_lazies = {}
        for (table, lazy_uuid_rows) in lazies.items():
            logger.debug("Making {} lazy proxies for objects in table '{}'"\
                         .format(len(lazy_uuid_rows), table))
            cls = self.serialization_schema.table_to_info[table].cls
            for row in lazy_uuid_rows:
                all_lazies[row.uuid] = self.make_lazy(cls, row.uuid)
        return all_lazies
