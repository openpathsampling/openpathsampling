import importlib
import re
import json
import networkx as nx
import numpy as np
import networkx.algorithms.dag as nx_dag
import collections
from tools import flatten_all, nested_update, group_by_function
from tools import is_iterable, is_mappable, is_numpy_iterable
import tools

# UUID recognition and encoding #####################################
# Things in here might be modified for performance optimization. In
# particular, it might be worth using a string representation of the UUID
# whenever possible (dicts with string keys have a special fast-path)


def has_uuid(obj):
    return hasattr(obj, '__uuid__')


def get_uuid(obj):
    # TODO: I can come up with a better string encoding than this
    try:
        return str(obj.__uuid__)
    except AttributeError as e:
        if obj is None:
            return obj
        else:
            raise e

def set_uuid(obj, uuid):
    obj.__uuid__ = long(uuid)


def encode_uuid(uuid):
    return "UUID(" + str(uuid) + ")"


def decode_uuid(uuid_str):
    return uuid_str[5:-1]


# use the regular expression when looking through an entire JSON string; use
# the is_uuid_string method for individual objects
encoded_uuid_re = re.compile("UUID\((?P<uuid>[0-9]+)\)")


def is_uuid_string(obj):
    return (
        isinstance(obj, (str, unicode))
        and obj[:5] == 'UUID(' and obj[-1] == ')'
    )


# Getting the list of UUIDs bsed on initial objets ###################
def caches_contain(key, cache_list):
    for cache in cache_list:
        if key in cache:
            return True
    return False

def filter_known_uuids(uuid_dict, cache_list):
    """Filters out UUIDs in the cache_list, returning what isn't cached"""
    return {uuid: value for (uuid, value) in uuid_dict.items()
            if not caches_contain(uuid, cache_list)}

class DefaultFindUUIDs(object):
    def __init__(self):
        self._mappable_classes = set()
        self._non_mappable_classes = set()
        self._iterable_classes = set()
        self._non_iterable_classes = set()

    @staticmethod
    def _get_with_cache(obj, key, true_set, false_set, check_method):
        # TODO: this is probably best moved to utils; I think we may use it
        # elsewhere as well
        if key in false_set:
            return False
        elif key in true_set:
            return True
        else:
            result = check_method(obj)
            set_for_obj = {True: true_set, False: false_set}[result]
            set_for_obj.add(key)
            return result

    def is_mappable(self, obj):
        return self._get_with_cache(
            obj=obj,
            key=obj.__class__,
            true_set=self._mappable_classes,
            false_set=self._non_mappable_classes,
            check_method=is_mappable
        )

    def is_iterable(self, obj):
        return self._get_with_cache(
            obj=obj,
            key=obj.__class__,
            true_set=self._iterable_classes,
            false_set=self._non_iterable_classes,
            check_method=lambda obj: \
                is_iterable(obj) and not is_numpy_iterable(obj)
        )

    # TODO: FIXME: add default_find_uuids as __call__ here; this should
    # really cut the time spent in is_iterable/is_mappable

def default_find_uuids(obj, cache_list):
    uuids = {}
    new_objects = []
    obj_uuid = get_uuid(obj) if has_uuid(obj) else None
    # filter known uuids: skip processing if known
    if caches_contain(obj_uuid, [uuids] + cache_list):
        return uuids, new_objects
    # UUID objects
    if obj_uuid:
        # print repr(obj)
        # print obj.to_dict().keys()
        uuids.update({obj_uuid: obj})
        new_objects.extend(list(obj.to_dict().values()))

    # mappables and iterables
    if is_mappable(obj):
        new_objects.extend([o for o in obj.keys() if has_uuid(o)])
        new_objects.extend(list(obj.values()))
    elif is_iterable(obj) and not is_numpy_iterable(obj):
        new_objects.extend(list(obj))
    return uuids, new_objects

# NOTE: this needs find everything, including if the iterable/mapping has a
# UUID, find that and things under it
def get_all_uuids(initial_object, known_uuids=None, class_info=None):
    """Find all UUID objects (to be stored)

    This searches through an initial object, finding *all* nested objects
    (including those in lists and dictionaries) that have UUIDs.

    Parameters
    ----------
    initial_object : object with UUID
        the object to search within
    known_uuids : dict of {uuid: object}
        objects that can be excluded from the search tree, presumably
        because they have already been searched and any object beneath them
        in the search tree also also already known
    class_info : :class:`.SerializationSchema`

    Returns
    -------
    dict of {uuid: object}
        objects found in the search
    """
    known_uuids = tools.none_to_default(known_uuids, {})
    objects = [initial_object]
    uuids = {}
    while objects:
        new_objects = []
        for obj in objects:
            # TODO: find a way to ensure that objects doesn't go over
            # duplicates here; see lprofile of default_find_uuids to see how
            # often abort due to being in cache in there, and whether we
            # should move the skip in here instead (how expensive is
            # info_from_instance?)
            info = class_info.info_from_instance(obj) \
                    if class_info else None
            if info:
                find_uuids = info.find_uuids
            else:
                find_uuids = default_find_uuids

            new_uuids, new_objs = find_uuids(obj=obj,
                                             cache_list=[uuids, known_uuids])

            uuids.update(new_uuids)
            new_objects.extend(new_objs)

        objects = new_objects
    return uuids


class SchemaFindUUIDs(object):
    def __init__(self, schema_entries):
        self.schema_entries = [
            (attr, attr_type) for (attr, attr_type) in schema_entries
            if attr_type in ['uuid', 'lazy', 'list_uuid']
        ]

    def __call__(self, obj, cache_list):
        uuids = {get_uuid(obj): obj}
        new_objects = []

        for (attr, attr_type) in self.schema_entries:
            attr_obj = getattr(obj, attr)
            if attr_type in ['uuid', 'lazy']:
                new_objects.append(attr_obj)
            elif attr_type == 'list_uuid':
                new_objects.extend(attr_obj)

        return uuids, new_objects


# TODO: I think this can be removed (only used by get_all_uuid_string)
def find_dependent_uuids(json_dct):
    dct = json.loads(json_dct)
    uuids = [decode_uuid(obj) for obj in flatten_all(dct)
             if is_uuid_string(obj)]
    return uuids

# TODO: I think this can be removed (not used?)
def get_all_uuid_strings(dct):
    all_uuids = []
    for uuid in dct:
        all_uuids.append(uuid)
        all_uuids += find_dependent_uuids(dct[uuid])
    return all_uuids


# NOTE: this only need to find until the first UUID: iterables/mapping with
# UUIDs aren't necessary here
def replace_uuid(obj, uuid_encoding):
    """Return storage-ready replacements for values in dict representation

    This is used by first creating the dict representation of an object by
    calling ``obj.to_dict()``. The resulting dict can be the ``obj``
    parameter of ``replace_uuid``. This algorithm is implemented
    recursively, so nested structures will call it again to create the
    correct nested locations of any UUID objects.

    Note that this does not go into the internal structure of any objects
    with UUIDs (such as UUID objects that are also iterable or mappable).
    The purpose here is to transform to the dict representation with UUIDs
    in the place of objects.

    Parameters
    ----------
    obj : object
        input; replace with UUID if is has one, or search for nested
        structures
    uuid_encoding : callable
        function that maps a UUID to a version that can be identified from
        the (JSON) serialized string representation of the object.

    Returns
    -------
    object
        input object with UUID objects replaced by the encoded form,
        leaving structure (of dicts, lists, etc) and all other objects the
        same
    """
    # this is UUID => string
    replacement = obj
    # fast exit for string keys
    if tools.is_string(obj):
        return replacement
    if has_uuid(obj):
        replacement = uuid_encoding(get_uuid(obj))
    elif is_mappable(obj):
        replacement = {
            replace_uuid(k, uuid_encoding): replace_uuid(v, uuid_encoding)
            for (k, v) in replacement.items()
        }
    elif is_iterable(obj) and not is_numpy_iterable(obj):
        replace_type = type(obj)
        replacement = replace_type([replace_uuid(o, uuid_encoding)
                                    for o in obj])
    return replacement


def to_dict_with_uuids(obj):
    dct = obj.to_dict()
    return replace_uuid(dct, uuid_encoding=encode_uuid)


def to_bare_json(obj):
    replaced = replace_uuid(obj, uuid_encoding=encode_uuid)
    return json.dumps(replaced)


def to_json_obj(obj):
    dct = to_dict_with_uuids(obj)
    dct.update({'__module__': obj.__class__.__module__,
                '__class__': obj.__class__.__name__})
    return json.dumps(dct)


def import_class(mod, cls):
    # TODO: this needs some error-checking
    mod = importlib.import_module(mod)
    cls = getattr(mod, cls)
    return cls

def search_caches(key, cache_list, raise_error=True):
    """Find UUID if it is in the cache_list dicts

    Parameters
    ----------
    key : str
        the UUID we're looking for
    cache_list : mapping or list of mapping
        caches that the objects are stored in (will be searched in order of
        the list). Mapping is {uuid: object}
    raise_error : bool
        whether to raise a KeyError if UUID not found; default True. If
        False, object not found returns None

    Returns
    -------
    object or None
        the object with the given UUID, or ``None`` if the object is not
        found and ``raise_error`` is ``False``.
    """
    if not isinstance(cache_list, list):
        cache_list = [cache_list]
    obj = None
    for cache in cache_list:
        if key in cache:
            obj = cache[key]
            break
    if obj is None and raise_error:
        raise KeyError("Missing key: " + str(key))
    return obj


def from_dict_with_uuids(obj, cache_list):
    """Replace encoded UUIDs with the actual objects.

    This is used to replace UUIDs as encoding for storage with actual
    objects, to complete reconstruction of the dict representation.  This
    method can be seen as the inverse process of :meth:`.replace_uuid`, an
    must be done before using the ``cls.from_dict()`` to properly
    instantiate the object.

    All input objects for the object being reconstructed must already be in
    the ``cache_list``. This means that the order is very important; that is
    controlled by :meth:`.get_reload_order`.

    Parameters
    ----------
    obj : object
        object to reconstruct; replace UUID strings with the actual objects
    cache_list : mapping or list of mapping
        existing objects, keyed by their UUIDs

    Returns
    -------
    object
        input object with UUID strings replaced by the actual objects
    """
    replacement = obj
    if is_uuid_string(obj):
        # raises KeyError if object hasn't been visited
        # (indicates problem in DAG reconstruction)
        uuid = decode_uuid(obj)
        replacement = search_caches(uuid, cache_list)
    elif tools.is_string(obj):
        # fast exit for string keys
        return obj
    elif is_mappable(obj):
        replacement = {from_dict_with_uuids(k, cache_list): \
                       from_dict_with_uuids(v, cache_list)
                       for (k, v) in obj.items()}
    elif is_iterable(obj) and not is_numpy_iterable(obj):
        replace_type = type(obj)
        replacement = replace_type([from_dict_with_uuids(o, cache_list)
                                    for o in obj])
    return replacement


def from_json_obj(uuid, table_row, cache_list):
    # NOTE: from_json only works with existing_uuids (DAG-ordering)
    dct = json.loads(table_row['json'])
    cls = import_class(dct.pop('__module__'), dct.pop('__class__'))
    dct = from_dict_with_uuids(dct, cache_list)
    obj = cls.from_dict(dct)
    set_uuid(obj, uuid)
    return obj


def _uuids_from_table_row(table_row, schema_entries):
    """Gather UUIDs from a table row (as provided by storage).

    This organizes the UUIDS that are included in the table row based on
    information from that row. It separated objects to be proxied ('lazy')
    from objects to be directly loaded ('uuid', 'list_uuid', 'json_obj').
    It also create the dependency dictionary (the entry for this row) that
    will be used to create the reconstruction DAG.

    This is for internal use; not to be part of the API.

    Parameters
    ----------
    table_row : object
        must have attributes as defined by ``schema_entries``, plus a
        ``uuid`` attribute. Typically comes directly from the backend.
    schema_entries : list of 2-tuple
        the pairs of (attribute_name, attribute_type) describing the columns
        from the ``table_row``. Should match the schema entry for the table
        that the table row comes from.

    Returns
    -------
    uuid : list
        list of UUIDs to be fully loaded
    lazy : set
        set of UUIDs to a lazy-loaded (i.e., as proxy)
    dependencies : dict
        length 1 dict mapping the input row's UUID to all UUIDs that it
        directly depends on (i.e., everything from ``uuid`` and ``lazy``).
    """
    # take the schema entries here, not the whole schema
    lazy = set([])
    uuid = set([])
    for (attr, attr_type) in schema_entries:
        if attr_type == 'uuid':
            uuid.add(getattr(table_row, attr))
        elif attr_type == 'list_uuid':
            # TODO: better to use encoded_uuid_re here?
            uuid_list = json.loads(getattr(table_row, attr))
            if uuid_list:
                # skip if None (or empty)
                uuid_list = {decode_uuid(u) for u in uuid_list}
                uuid.update(uuid_list)
        elif attr_type == 'json_obj':
            json_dct = getattr(table_row, attr)
            uuid.update(set(encoded_uuid_re.findall(json_dct)))
        elif attr_type == 'lazy':
            lazy.add(getattr(table_row, attr))
        # other cases aren't UUIDs and are ignored
    # remove all cases of None as a UUID to depend on
    # TODO: should None be in the UUID list even?
    # TODO: can we return the set here?
    dependencies = {table_row.uuid: (uuid | lazy) - {None}}
    return (list(uuid), lazy, dependencies)


def get_all_uuids_loading(uuid_list, backend, schema, existing_uuids=None):
    """Get all information to reload from UUIDs.

    This is the main function for identifying objects to reload from
    storage. It returns the table rows to load (sorted by table), the UUIDs
    of objects to lazy-load (sorted by table), and the dictionary of
    dependencies, which can be used to create the reconstruction DAG.
    """
    if existing_uuids is None:
        existing_uuids = {}
    known_uuids = set(existing_uuids.keys())
    uuid_to_table = {}
    all_table_rows = []
    lazy = set([])
    dependencies = {}
    while uuid_list:
        new_uuids = {uuid for uuid in uuid_list if uuid not in known_uuids}
        uuid_rows = backend.load_uuids_table(new_uuids)
        new_table_rows = backend.load_table_data(uuid_rows)
        uuid_to_table.update({r.uuid: backend.uuid_row_to_table_name(r)
                              for r in uuid_rows})

        uuid_list = []
        for row in new_table_rows:
            entries = schema[uuid_to_table[row.uuid]]
            loc_uuid, loc_lazy, deps = _uuids_from_table_row(row, entries)
            uuid_list += loc_uuid
            lazy.update(loc_lazy)
            dependencies.update(deps)

        all_table_rows += new_table_rows
        known_uuids |= new_uuids
        uuid_list = {uuid for uuid in uuid_list if uuid not in known_uuids}

    return (all_table_rows, lazy, dependencies, uuid_to_table)


def dependency_dag(dependent_uuids, dag=None):
    if dag is None:
        dag = nx.DiGraph()
    for from_node, to_nodes in dependent_uuids.items():
        if to_nodes:
            dag.add_edges_from([(from_node, to_node)
                                for to_node in to_nodes])
    if not nx_dag.is_directed_acyclic_graph(dag):
        raise RuntimeError("Reconstruction DAG not acyclic?!?!")
    return dag

def dag_reload_order(dag):
    return list(reversed(list(nx_dag.topological_sort(dag))))

def get_reload_order(to_load, dependencies):
    dag = dependency_dag(dependencies)
    no_deps = {row.uuid for row in to_load}
    no_deps.difference_update(set(dag.nodes))
    ordered_uuids = list(no_deps) + dag_reload_order(dag)
    return ordered_uuids
