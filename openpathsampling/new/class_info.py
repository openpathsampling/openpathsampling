import tools

from serialization import DefaultSerializer, DefaultDeserializer

import cloudpickle

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
                 lookup_result=None):
        self.table = table
        self.cls = cls
        self.serializer = serializer
        self.deserializer = deserializer
        if lookup_result is None:
            lookup_result = cls
        self.lookup_result = lookup_result

    def set_defaults(self, schema):
        if self.serializer is None:
            self.serializer = DefaultSerializer(schema, self.table,
                                                self.cls)
        if self.deserializer is None:
            self.deserializer = DefaultDeserializer(schema, self.table,
                                                    self.cls)

    def __repr__(self):
        return ("ClassInfo(table=" + self.table + ", cls=" + str(self.cls)
                + ", lookup_result=" + str(self.lookup_result) + ")")


class ClassInfoContainer(object):
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
    def __init__(self, default_info, class_info_list=None):
        class_info_list = tools.none_to_default(class_info_list, [])
        self.lookup_to_info = {}
        self.table_to_info = {}
        self.class_info_list = []
        self.default_info = default_info
        self.missing_table = ClassInfo(table="__missing__", cls=None)
        self.add_class_info(default_info)
        for info in class_info_list:
            self.add_class_info(info)

    def add_class_info(self, info_node):
        # check that we're not in any existing
        self.class_info_list.append(info_node)
        # TODO: is it possible to generalize the special cases (use the type
        # of the expected return to identify the function to call -- won't
        # be completely general, but may simplify some)
        self.lookup_to_info.update({info_node.lookup_result: info_node})
        self.table_to_info.update({info_node.table: info_node})

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

    def __repr__(self):  # pragma: no cover
        return ("ClassInfoContainer(default_info=" + repr(self.default_info)
                + ", class_info_list=" + repr(self.class_info_list) + ")")


