import os
import collections
from collections import abc

from .storage import universal_schema

import logging
logger = logging.getLogger(__name__)

UUIDTableRow = collections.namedtuple("UUIDTableRow", ['uuid', 'table_name'])

class MemoryStorageTable(abc.MutableMapping):
    """Simple wrapper for a Dict; transforms Dict values to namedtuple
    """
    def __init__(self, table_name, schema_entries):
        self.table_name = table_name
        self.schema_entries = schema_entries
        self.Namespace = self._make_namespace(table_name, schema_entries)
        self._dict = {}

    def __reduce__(self):
        # remove dynamically created namedtuples
        return (
            self.__class__._from_reduce_tuple, (
                self.table_name,
                self.schema_entries,
                {key: tuple(val) for key, val in self._dict.items()}
            )
        )

    @classmethod
    def _from_reduce_tuple(cls, table_name, schema_entries, dct):
        obj = cls(table_name, schema_entries)
        fields = ['uuid'] + [name for name, _ in schema_entries]
        for key, val in dct.items():
            obj[key] = {key: v for key, v in zip(fields, val)}
        return obj

    @staticmethod
    def _make_namespace(name, schema_entries):
        attrs = [attr for attr, _ in schema_entries]
        namespace = collections.namedtuple(name, ['uuid'] + attrs)
        return namespace

    def __getitem__(self, key):
        return self._dict[key]

    def __setitem__(self, key, value):
        self._dict[key] = self.Namespace(**value)

    def __delitem__(self, key):
        del self._dict[key]

    def __iter__(self):
        return iter(self._dict)

    def __len__(self):
        return len(self._dict)

    def copy(self):
        my_copy = MemoryStorageTable(self.table_name, self.schema_entries)
        my_copy._dict = self._dict.copy()
        return my_copy


class MemoryStorageBackend(object):
    """Storage backend for memory.

    This is primarily used for testing and for serialization over a network.
    This can be considered as a simplest working example of a backend, and
    therefore may be useful for developers.
    """
    def __init__(self, mode='w'):
        self.mode = mode  # technically always 'a', but setup may differ
        self.filename = "<Memory>"
        self.schema = {}
        self.namedtuples = {}  # maps table name to associated namedtuple
        self.data = {}  # maps table name to Dict[UUID, namedtuple]
        self.uuid_table = {}  # maps UUID to table name
        self._table_to_class = {}

    @property
    def identifier(self):
        return hex(id(self))

    def close(self):
        pass  # no such thing as closing here, but may be needed for API

    def register_type(self, type_str, backend_type):
        pass  # no need to do anything here?

    def register_schema(self, schema, table_to_class, metadata=None):
        for table_name, schema_entries in schema.items():
            logger.debug("Registering table: %s" % table_name)
            self.data[table_name] = MemoryStorageTable(table_name,
                                                       schema_entries)

        self.schema.update(schema)
        self._table_to_class.update(table_to_class)

    def has_table(self, table_name):
        return table_name in self.data

    def register_storable_function(self, table_name, result_type):
        self.data[table_name] = {}

    def add_storable_function_results(self, table_name, result_dict):
        # TODO: safety checks that we don't overwrite existing
        self.data[table_name].update(result_dict)

    def load_storable_function_results(self, table_name, uuids):
        return {uuid: self.data[table_name][uuid] for uuid in uuids}

    def load_storable_function_table(self, table_name):
        return self.data[table_name].copy()

    def add_to_table(self, table_name, objects):
        uuids = [obj['uuid'] for obj in objects]
        # TODO: check for existence (for safety)
        self.uuid_table.update({uuid: table_name for uuid in uuids})
        self.data[table_name].update({obj['uuid']: obj for obj in objects})

    def load_n_rows_from_table(self, table_name, first_row, n_rows):
        # not the prettiest implementation, but everything must already be
        # in memory here
        rows = list(self.data[table_name].values())
        return rows[first_row:first_row + n_rows]

    def load_uuids_table(self, uuids, ignore_missing=False):
        exists = [uuid for uuid in uuids if uuid in self.uuid_table]
        not_exists = set(uuids) - set(exists)
        if not ignore_missing and len(not_exists) > 0:
            # at some point we need to raise errors here
            pass

        results = [UUIDTableRow(uuid, self.uuid_table[uuid])
                   for uuid in exists]
        return results

    def load_table_data(self, uuid_table_rows):
        results = [self.data[row.table_name][row.uuid]
                   for row in uuid_table_rows]
        return results

    def database_schema(self):
        return self.schema

    def get_representative(self, table_name):
        return next(iter(self.data[table_name].values()))

    @property
    def table_to_class(self):
        return self._table_to_class

    def uuid_row_to_table_name(self, uuid_row):
        return uuid_row.table_name

    def table_iterator(self, table_name):
        return iter(self.data[table_name].values())

    def table_len(self, table_name):
        return len(self.data[table_name])

    def table_get_item(self, table_name, item):
        return list(self.data[table_name].values())[item]
