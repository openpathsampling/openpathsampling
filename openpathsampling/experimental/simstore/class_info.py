import logging

from . import tools
from .serialization import SchemaSerializer, SchemaDeserializer
from .serialization_helpers import SchemaFindUUIDs, has_uuid
from .serialization_helpers import encoded_uuid_re, get_reload_order
from .serialization_helpers import get_all_uuids
from .my_types import uuid_types, uuid_list_types, json_obj_types

from . import attribute_handlers

import json


logger = logging.getLogger(__name__)


class ClassInfo(object):
    """
    Parameters
    ----------
    table : str
        the table name for this object type
    cls : class
        the class for this object type (used in the default deserializer)
    serializer : callable
        serializer for this object type (None gives the default serializer)
    deserializer : callable
        deserializer for this object type (None gives the default
        deserializer)
    lookup_result : any
        the result when ClassInfoContainer looks up objects in this table
    find_uuids : callable
        a shortcut for finding objects with UUIDs contained within this
        object
    """
    def __init__(self, table, cls, serializer=None, deserializer=None,
                 lookup_result=None, find_uuids=None,
                 safe_deserializer=None):
        self.table = table
        self.cls = cls
        self.serializer = serializer
        self.deserializer = None
        self.unsafe_deserializer = deserializer
        self.safe_deserializer = safe_deserializer
        if lookup_result is None:
            lookup_result = cls
        self.lookup_result = lookup_result
        self.find_uuids = find_uuids

    def set_defaults(self, schema, handlers):
        table = self.table if self.table in schema else None
        self.serializer = tools.none_to_default(
            self.serializer,
            SchemaSerializer(schema, table, self.cls, handlers)
        )
        self.unsafe_deserializer = tools.none_to_default(
            self.unsafe_deserializer,
            SchemaDeserializer(schema, table, self.cls, handlers)
        )
        self.safe_deserializer = tools.none_to_default(
            self.safe_deserializer,
            self.unsafe_deserializer
        )
        default_entries = schema.get(self.table, [])
        self.find_uuids = tools.none_to_default(
            self.find_uuids, SchemaFindUUIDs(default_entries)
        )
        self.set_safemode(True)

    def set_safemode(self, mode):
        self.deserializer = {True: self.safe_deserializer,
                             False: self.unsafe_deserializer}[mode]

    def __repr__(self):  # pragma: no cover
        return ("ClassInfo(table=" + self.table + ", cls=" + str(self.cls)
                + ", lookup_result=" + str(self.lookup_result)
                + ", find_uuids=" + str(self.find_uuids) + ")")


class SerializationSchema(object):
    """Connect table to serialization information.

    This is the primary location for information connecting application
    objects (i.e., OPS) to the serialization mechanisms. In particular, the
    __getattr__ in this can take either an instance or a table name, and in
    either case links it to the ClassInfo. This lets us link an object
    instance to the name of the table where it should be stored, and to
    link the table to the serialization methods.

    Notes
    -----
    Linking instances to tables can require special treatment. The default
    behavior is to use the class of the object to identify the table, and
    this works for most cases. However, sometimes a class can map to more
    than one table, and therefore is not a unique key. For example, a single
    class might be used to represent data with different dimensions, and
    therefore require different tables (e.g, coordinates for different
    systems). In such cases, the SerializationSchema needs to be subclassed
    with specialized information.
    """
    def __init__(self, default_info, sfr_info=None, schema=None,
                 class_info_list=None, handlers=None):
        class_info_list = tools.none_to_default(class_info_list, [])
        handlers = tools.none_to_default(handlers,
                                         attribute_handlers.DEFAULT_HANDLERS)
        self.attribute_handlers = handlers
        self.schema = {}
        self.lookup_to_info = {}
        self.table_to_info = {}
        self.class_info_list = []
        self.default_info = default_info
        self.sfr_info = sfr_info
        self.missing_table = ClassInfo(table="__missing__", cls=None)
        self.add_class_info(default_info)
        if sfr_info is not None:
            self.add_class_info(sfr_info)
        self.register_info(class_info_list, schema)

    # TODO: I think that this can be made private; used by __init__ and
    # within register_info
    def add_class_info(self, info_node):
        # check that we're not in any existing
        self.class_info_list.append(info_node)
        # TODO: is it possible to generalize the special cases (use the type
        # of the expected return to identify the function to call -- won't
        # be completely general, but may simplify some)
        self.lookup_to_info.update({info_node.lookup_result: info_node})
        immutable_tables = [self.default_info.table,
                            self.missing_table.table]
        sfr_table_as_list = (
            [self.sfr_info.table] if self.sfr_info is not None else []
        )
        reserved = info_node.table in immutable_tables + sfr_table_as_list
        exists = info_node.table in self.table_to_info
        if not (reserved and exists):
            self.table_to_info.update({info_node.table: info_node})

    def copy(self):
        # NOTE: subclasses must override if __init__ sig changes
        dup = self.__class__(
            default_info=self.default_info,
            sfr_info=self.sfr_info,
            schema=self.schema,
            class_info_list=self.class_info_list
        )
        return dup

    def backend_type(self, type_name):
        """Obtain the general type name from a specific type string

        Example: this tells that no matter what the specific shape of a
        numpy array is, the backend should register it as a NumPy array. In
        addition, this can (in the future will) provide a length value to
        tell the backend how many bytes to allocate for an array.

        Parameters
        ----------
        type_name: str
            specfic type string for an attribute

        Returns
        -------
        Tuple[str, Any]
            the type name as known to the backend, and the length
            information for fixed length columns (length aspect not yet
            implemented)
        """
        for handler in self.attribute_handlers:
            if handler.is_my_type(type_name):
                handler_obj = handler.from_type_string(type_name)
                return handler_obj.backend_type, handler_obj.type_size
        # current default is to return the input; may change
        return type_name, None

    def register_info(self, class_info_list, schema=None):
        schema = tools.none_to_default(schema, {})
        self.schema.update(schema)
        for info in class_info_list:
            info.set_defaults(schema, self.attribute_handlers)
            self.add_class_info(info)

    def set_safemode(self, mode):
        for class_info in self.class_info_list:
            class_info.set_safemode(mode)

    @property
    def tables(self):
        """list of tables from the included class info objects"""
        tables = [info.table for info in self.class_info_list]
        tables.append(self.default_info.table)
        return tables

    def lookup_key(self, item):
        if not self.is_special(item):
            lookup_key = item.__class__
        else:
            lookup_key = self.special_lookup_key(item)
        return lookup_key

    def is_special(self, item):
        return False

    def special_lookup_key(self, item):
        raise NotImplementedError("No special types implemented")

    def get_special(self, item):
        lookup = self.special_lookup_key(item)
        try:
            return self.lookup_to_info[lookup]
        except KeyError:
            return self.missing_table

    def __getitem__(self, item):
        if tools.is_string(item):
            return self.table_to_info[item]
        else:
            info = self.info_from_instance(item)
            if info is None:
                raise KeyError("'%s'" % repr(item))
            else:
                return info

        # elif self.is_special(item):
            # return self.get_special(item)
        # else:
            # lookup = self.lookup_key(item)
            # try:
                # return self.lookup_to_info[lookup]
            # except KeyError as e:
                # if isinstance(item, self.default_info.cls):
                    # return self.default_info
                # else:
                    # raise e

    def info_from_instance(self, item):
        if not has_uuid(item):
            return None
        if self.is_special(item):
            return self.get_special(item)
        else:
            lookup = self.lookup_key(item)
            if lookup in self.lookup_to_info:
                return self.lookup_to_info[lookup]
            elif isinstance(item, self.default_info.cls):
                return self.default_info
            else:
                return None

    def add_missing_table_from_instance(self, lookup, obj):
        raise NotImplementedError("No special types implemented")

    # TODO: this is currently done in storage; should it be here?
    def _missing_table_update(self, by_table):
        missing = by_table.pop('__missing__')
        logger.info("Identifying tables for %d objects of unknown "
                     + "table type", len(missing))

        lookup_examples = set([])
        for obj in missing.values():
            lookup = self.lookup_key(obj)
            if lookup not in lookup_examples:
                self.add_missing_table_from_instance(lookup, obj)
                lookup_examples |= {lookup}

        get_table_name = lambda uuid, obj_: self[obj_].table
        missing_by_table = tools.dict_group_by(missing, get_table_name)
        logger.info("Identified %d new tables: %s", len(missing_by_table),
                    str(list(missing_by_table.keys())))
        by_table.update(missing_by_table)
        return by_table

    # def serialize(self, obj, storage=None):
        # """
        # Serialize all objects in ``obj`` to table-grouped JSON-ready dict.

        # Notes
        # -----

        # This can also take multiple objects, but they must be collected in a
        # hashable container (e.g., tuple, not a list).

        # Parameters
        # ----------
        # obj : object
            # the object to be serialized
        # storage : :class:`.Storage`
            # storage object with connection to backend and cache (if
            # relevant)

        # Returns
        # -------
        # dict :
            # {table: {uuid: {attr: value}}}
        # """
        # # TODO: correct the return type (no dict of uuid to description,
        # # just a list of dicts for each object, including uuid as attr
        # logger.debug("Starting serialization")
        # results = {}
        # cache = {} if not storage else storage.cache
        # if storage and storage.uuids_in_storage([get_uuid(obj)]):
            # return  # return what? is None right?

        # logger.debug("Listing all included objects to serialize")
        # uuids = get_all_uuids(obj, cache)  # TODO: replace w SchemaFindUUIDs

        # if storage:
            # exists = storage.uuids_in_storage(list(uuids.keys()))
            # for existing in exists:
                # del uuids[existing.uuid]

        # get_table_name = lambda uuid, obj_: self[obj_].table
        # by_table = tools.dict_group_by(uuids, key_extract=get_table_name)

        # # fix for missing tables
        # if '__missing__' in by_table:
            # by_table = self._missing_table_update(by_table)
            # # NOTE: if using a storage, this will still need to register the
            # # tables with the backend -- that should be handled as part of
            # # the save process TODO

        # logger.debug("Serializing objects from %d tables: %s",
                     # len(by_table), str(list(by_table.keys())))
        # serialized_by_table = {}
        # for (table, table_uuids) in by_table.items():
            # logger.debug("Serializing %d objects from table %s",
                         # len(by_table[table]), table)
            # serialize = self[table].serializer
            # serialized_table = []
            # for o in table_uuids.values():
                # if table == 'simulation_objects':
                    # logger.debug(str(o))
                # serialized_table.append(serialize(o))
            # serialized_by_table[table] = serialized_table
            # # serialized_by_table[table] = [serialize(o)
                                          # # for o in table_uuids.values()]

        # return serialized_by_table


    # def deserialize(self, serialized_by_table, known_uuids=None):
        # # 1. generate dependencies lists
        # # 2. get reload order
        # # 3. reconstruct UUIDs
        # known_uuids = tools.none_to_default(known_uuids, {})
        # dependencies = {}
        # to_load = {}
        # uuid_to_table = {}
        # for (table, object_list) in serialized_by_table.items():
            # schema_entries = self.schema[table]
            # logger.debug("Restoring %d objects from table %s",
                         # len(object_list), table)
            # for item_dct in object_list:
                # uuid = item_dct['uuid']
                # uuid_to_table[uuid] = table
                # to_load[uuid] = tools.SimpleNamespace(**item_dct)
                # if table == self.default_info.table:
                    # item_json = json.dumps(item_dct)
                    # dependencies[uuid] = set(encoded_uuid_re.findall(item_json))
                # else:
                    # # TODO: this can be cleaned up with a dict mapping types
                    # # to method for finding dependencies
                    # # find_deps is a defaultdict mapping to appropriate
                    # # functions, default is fcn returning empty set
                    # # for (entry, e_type) in schema_entries:
                        # # deps |= find_deps[e_type](entry)
                    # uuid_entries = [entry
                                    # for (entry, e_type) in schema_entries
                                    # if e_type in uuid_types]
                    # deps = set([item_dct[entry] for entry in uuid_entries])
                    # uuid_list_entries = [entry
                                         # for (entry, e_type) in schema_entries
                                         # if e_type in uuid_list_types]
                    # deps |= set(sum([encoded_uuid_re.findall(item_dct[entry])
                                     # for entry in uuid_list_entries], []))

                    # json_entries = [entry
                                    # for (entry, e_type) in schema_entries
                                    # if e_type in json_obj_types]
                    # for json_entry in json_entries:
                        # json_val = item_dct[json_entry]
                        # deps |= set(encoded_uuid_re.findall(json_val))

                    # dependencies[uuid] = deps

        # ordered_uuids = get_reload_order(list(to_load.values()), dependencies)
        # logger.debug("Restore order: %s", str(ordered_uuids))
        # return self.reconstruct_uuids(ordered_uuids, uuid_to_table, to_load,
                                      # known_uuids)


    # NOTE: this doesn't seem to be used yet.... should we move
    # functionality here or keep in storage.deserialize_uuids?
    # right now, leaning toward storage....
    def reconstruct_uuids(self, ordered_uuids, uuid_to_table, to_load,
                          known_uuids=None, new_uuids=None):
        """
        Parameters
        ----------
        ordered_uuids
        uuid_to_table
        known_uuids : dict {uuid: object}
            immutable list of known uuids
        new_uuids : dict {uuid: object}
            mutable list of uuids that have been reloaded, updated after
            these uuids have been reloaded (also returned)
        """
        # (this is just the current storage.deserialize_uuids)
        known_uuids = tools.none_to_default(known_uuids, {})
        new_uuids = tools.none_to_default(new_uuids, {})
        logger.debug("Reconstructing from %d objects", len(ordered_uuids))
        for uuid in ordered_uuids:
            if uuid not in known_uuids and uuid not in new_uuids:
                table = uuid_to_table[uuid]
                reconstruct = self[table].deserializer
                obj = reconstruct(uuid, to_load[uuid],
                                  [new_uuids, known_uuids])
                new_uuids[uuid] = obj

        return new_uuids


    def __repr__(self):  # pragma: no cover
        return ("SerializationSchema(default_info=" + repr(self.default_info)
                + ", class_info_list=" + repr(self.class_info_list) + ")")


# TODO: remove this alias
ClassInfoContainer = SerializationSchema
