from .wrapper import SimStoreWrapper
from .serialization_helpers import get_uuid

from collections import abc

class TaggedObject(SimStoreWrapper):
    pass  # subclass simply so it can be recognized distinct from wrapper


class TagsTable(abc.Mapping):
    """StorageTable-like object for tags.

    Tags allow storage of arbitrary objects. The ``tags`` table can contain
    any JSON-serializable object, or any UUID-containing object that can
    otherwise be stored by SimStore. Usage is dict-like, using a string name
    for the tag, e.g., ``tags['initial_conditions'] = data`` to set or
    ``data = tags['initial_conditions']`` to load.

    Parameters
    ----------
    storage: :class:`.GeneralStorage`
        the storage object associated with this table
    """
    # we inherit from Mapping, not MutableMapping, because delitem, pop,
    # popitem, etc don't make sense here
    def __init__(self, storage):
        self.storage = storage
        self.table = 'tags'
        self.tagged_objects = self._load_tags()

    def _load_tags(self):
        # TODO: maybe I should just make the table? Need better way for that
        if self.storage.backend.has_table(self.table):
            backend_iter = self.storage.backend.table_iterator(self.table)
            tagged_objects = {row.name: row.content for row in backend_iter}
        else:
            tagged_objects = {}
        return tagged_objects

    def __iter__(self):
        for name in self.tagged_objects:
            yield self[name]

    def __len__(self):
        return len(self.tagged_objects)

    def keys(self):
        return self.tagged_objects.keys()

    def __getitem__(self, tag):
        uuid = self.tagged_objects[tag]
        wrapped = self.storage.load([uuid])[0]
        return wrapped.content

    def __setitem__(self, tag, obj):
        if tag in self.tagged_objects:
            raise RuntimeError("A tag named '%s' already exists." % tag)
        wrapped = TaggedObject(obj).named(tag)
        self.storage.save(wrapped)
        wrapped_uuid = get_uuid(wrapped)
        self.tagged_objects[tag] = wrapped_uuid
        self.storage.backend.add_tag(self.table, tag, wrapped_uuid)
