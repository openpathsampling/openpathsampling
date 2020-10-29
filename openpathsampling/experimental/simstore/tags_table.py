from .wrapper import SimStoreWrapper
from .storage import StorageTable
from .serialization_helpers import get_uuid

from collections import abc

class TaggedObject(SimStoreWrapper):
    pass  # subclass simply so it can be recognized distinct from wrapper


class TagsTable(abc.Sequence):
    # we inherit from Sequence, not MutableSequence, because delitem and
    # insert don't make sense here
    def __init__(self, storage):
        self.storage = storage
        self.table = 'tags'
        self.tagged_objects = self._load_tags()

    def _load_tags(self):
        backend_iter = self.storage.backend.table_iterator(self.table)
        tagged_objects = {row.name: row.content for row in backend_iter}
        return tagged_objects

    def __len__(self):
        return len(self.tagged_objects)

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
