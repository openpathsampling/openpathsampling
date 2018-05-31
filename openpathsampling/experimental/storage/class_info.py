import tools

from serialization import DefaultSerializer, DefaultDeserializer
from serialization_helpers import SchemaFindUUIDs, has_uuid

try:
    import ujson as json
except ImportError:
    import json

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
    """
    def __init__(self, table, cls, serializer=None, deserializer=None,
                 lookup_result=None, find_uuids=None):
        self.table = table
        self.cls = cls
        self.serializer = serializer
        self.deserializer = deserializer
        if lookup_result is None:
            lookup_result = cls
        self.lookup_result = lookup_result
        self.find_uuids = find_uuids

    def set_defaults(self, schema):
        self.serializer = tools.none_to_default(
            self.serializer,
            DefaultSerializer(schema, self.table, self.cls)
        )
        self.deserializer = tools.none_to_default(
            self.deserializer,
            DefaultDeserializer(schema, self.table, self.cls)
        )
        self.find_uuids = tools.none_to_default(
            self.find_uuids,
            SchemaFindUUIDs(schema[self.table])
        )


    def __repr__(self):
        return ("ClassInfo(table=" + self.table + ", cls=" + str(self.cls)
                + ", lookup_result=" + str(self.lookup_result) + ")")



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
    systems). In such cases, the ClassInfoContainer needs to be subclassed
    with specialized information.
    """
    def __init__(self, default_info, schema=None, class_info_list=None):
        class_info_list = tools.none_to_default(class_info_list, [])
        self.schema = {}
        self.lookup_to_info = {}
        self.table_to_info = {}
        self.class_info_list = []
        self.default_info = default_info
        self.missing_table = ClassInfo(table="__missing__", cls=None)
        self.add_class_info(default_info)
        self.register_info(class_info_list, schema)

    def add_class_info(self, info_node):
        # check that we're not in any existing
        self.class_info_list.append(info_node)
        # TODO: is it possible to generalize the special cases (use the type
        # of the expected return to identify the function to call -- won't
        # be completely general, but may simplify some)
        self.lookup_to_info.update({info_node.lookup_result: info_node})
        self.table_to_info.update({info_node.table: info_node})

    def register_info(self, class_info_list, schema=None):
        schema = tools.none_to_default(schema, {})
        self.schema.update(schema)
        for info in class_info_list:
            info.set_defaults(schema)
            self.add_class_info(info)

    @property
    def tables(self):
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
        return NotImplementedError("No special types implemented")

    def get_special(self, item):
        lookup = self.special_lookup_key(item)
        try:
            return self.lookup_to_info[lookup]
        except KeyError:
            return self.missing_table

    def __getitem__(self, item):
        # TODO: base this off of info_from_instance
        if tools.is_string(item):
            return self.table_to_info[item]
        elif self.is_special(item):
            return self.get_special(item)
        else:
            lookup = self.lookup_key(item)
            try:
                return self.lookup_to_info[lookup]
            except KeyError as e:
                if isinstance(item, self.default_info.cls):
                    return self.default_info
                else:
                    raise e

    def info_from_instance(self, item):
        if not has_uuid(item):
            return None
        if self.is_special(item):
            self.get_special(item)
        else:
            lookup = self.lookup_key(item)
            if lookup in self.lookup_to_info:
                return self.lookup_to_info[lookup]
            elif isinstance(item, self.default_info.cls):
                return self.default_info
            else:
                return None


    def serialize(self, obj, storage=None):
        """
        Serialize all objects in ``obj`` to table-grouped JSON-ready dict.

        Parameters
        ----------
        obj : object
            the object to be serialized
        storage : :class:`.Storage`
            storage object with connection to backend and cache (if
            relevant)

        Returns
        -------
        dict :
            {table: {uuid: {attr: value}}}
        """
        logger.debug("Starting serialization")
        results = {}
        cache = {} if not storage else storage.cache
        if storage and storage.uuids_in_storage([get_uuid(obj)]):
            return  # return what? is None right?

        logger.debug("Listing all included objects to serialize")
        uuids = get_all_uuids(obj, cache)  # TODO: replace this

        if storage:
            exists = storage.uuids_in_storage(list(uuids.keys()))
            for existing in exists:
                del uuids[existing.uuid]

        get_table_name = lambda uuid, obj_: self[obj_].table
        by_table = tools.dict_group_by(uuids, key_extract=get_table_name)

        # fix for missing tables

        logger.debug("Serializing objects from %d tables: %s",
                     len(by_table), str(list(by_table.keys())))
        serialized_by_table = {}
        for (table, objects) in by_table.items():
            logger_debug("Serializing %d objects from table %s",
                         len(by_table[table]), table)
            serialize = self[table].serializer
            serialized_by_table[table] = [serialize(o) for o in objects]

        return serialized_by_table


    def deserialize(self, serialized_by_table):
        # 1. generate dependencies lists
        # 2. get reload order
        # 3. reconstruct UUIDs
        dependencies = {}
        to_load = []
        uuid_to_table = {}
        for (table, objects) in serialized_by_table.items():
            schema_entries = self.schema[table]
            for (uuid, item_dct) in objects.items():
                uuid_to_table[uuid] = table
                to_load.append(tools.SimpleNamespace(**item_dct))
                item_json = json.dumps(item_dct)
                dependencies[uuid] = set(encoded_uuid_re.findall(item_json))

        # TODO: this is the reload order now
        ordered_uuids = get_reload_order(to_load, dependencies)
        # TODO: self.reconstruct_uuids(ordered_uuids, uuid_to_table)
        # (this is just the current storage.deserialize_uuids)




    def __repr__(self):  # pragma: no cover
        return ("SerializationSchema(default_info=" + repr(self.default_info)
                + ", class_info_list=" + repr(self.class_info_list) + ")")


ClassInfoContainer = SerializationSchema
