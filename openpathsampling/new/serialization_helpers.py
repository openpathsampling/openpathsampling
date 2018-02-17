import importlib
import re
import ujson
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


encoded_uuid_re = re.compile("UUID\((?P<uuid>[0-9]+)\)")


def is_uuid_string(obj):
    return (
        isinstance(obj, (str, unicode))
        and obj[:5] == 'UUID(' and obj[-1] == ')'
    )


# Getting the list of UUIDs bsed on initial objets ###################

# NOTE: this needs find everything, including if the iterable/mapping has a
# UUID, find that and things under it
def get_all_uuids(initial_object, known_uuids=None):
    known_uuids = tools.none_to_default(known_uuids, {})
    objects = {initial_object}
    uuids = {}
    uuid_or_none = lambda o: get_uuid(o) if has_uuid(o) else None
    while objects:
        new_objects = []
        for obj in objects:
            obj_uuid = uuid_or_none(obj)
            # filter known uuids: skip processing if known
            if obj_uuid in uuids or obj_uuid in known_uuids:
                continue
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
        objects = new_objects
    return uuids


# TODO: I think this can be removed (only used by get_all_uuid_string)
def find_dependent_uuids(json_dct):
    dct = ujson.loads(json_dct)
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
    return ujson.dumps(replaced)


def to_json_obj(obj):
    dct = to_dict_with_uuids(obj)
    dct.update({'__module__': obj.__class__.__module__,
                '__class__': obj.__class__.__name__})
    return ujson.dumps(dct)


def import_class(mod, cls):
    # TODO: this needs some error-checking
    mod = importlib.import_module(mod)
    cls = getattr(mod, cls)
    return cls

def search_caches(key, cache_list, raise_error=True):
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
    dct = ujson.loads(table_row['json'])
    cls = import_class(dct.pop('__module__'), dct.pop('__class__'))
    dct = from_dict_with_uuids(dct, cache_list)
    obj = cls.from_dict(dct)
    set_uuid(obj, uuid)
    return obj


def uuids_from_table_row(table_row, schema_entries):
    # take the schema entries here, not the whole schema
    lazy = set([])
    uuid = []
    for (attr, attr_type) in schema_entries:
        if attr_type == 'uuid':
            uuid.append(getattr(table_row, attr))
        elif attr_type == 'list_uuid':
            # TODO: can find_dependent_uuids work here?
            uuid_list = ujson.loads(getattr(table_row, attr))
            if uuid_list:
                # skip if None (or empty)
                uuid_list = [decode_uuid(u) for u in uuid_list]
                uuid.extend(uuid_list)
        elif attr_type == 'json_obj':
            json_dct = getattr(table_row, attr)
            uuid.extend(encoded_uuid_re.findall(json_dct))
        elif attr_type == 'lazy':
            lazy.add(getattr(table_row, attr))
        # other cases aren't UUIDs and are ignored
    # remove all existiences of None as a UUID to depend on
    dependencies = {table_row.uuid: (set(uuid) | lazy) - {None}}
    return (uuid, lazy, dependencies)


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
            loc_uuid, loc_lazy, deps = uuids_from_table_row(row, entries)
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
