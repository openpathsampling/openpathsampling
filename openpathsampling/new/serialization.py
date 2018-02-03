import importlib
import collections
import ujson
import networkx as nx
import networkx.algorithms.dag as nx_dag
from tools import flatten_all, nested_update

# mappable/iterable identification ##################################
# TODO: may just hard-code these; this seems to be the proper way

def is_mappable(obj):
    return isinstance(obj, collections.Mapping)

def is_iterable(obj):
    return isinstance(obj, collections.Iterable)

# UUID recognition and encoding #####################################

def uuid_test(obj):
    # TODO: or isinstance? or other? try a few for performance checks
    return hasattr(obj, '__uuid__')

# TODO: try a few UUID encodings for performance
def encode_uuid(uuid):
    return "UUID(" + str(uuid) + ")"

def decode_uuid(uuid_str):
    return long(uuid_str[5:-1])

def is_uuid_string(obj):
    return (
        isinstance(obj, (str, unicode))
        and obj[:5] == 'UUID(' and obj[-1] == ')'
    )

# TODO: have a special UUID encoding for dict keys? string keys have special
# fast-path


# Getting the list of UUIDs bsed on initial objets ###################

# TODO: does this work with arbitrary nested yet? (flatten_all?)
# NOTE: this needs find everything, including if the iterable/mapping has a
# UUID, find that and things under it
def get_all_uuids(initial_object):
    uuid_dict = {initial_object.__uuid__: initial_object}
    with_uuid = [o for o in flatten_all(initial_object.to_dict())
                 if uuid_test(o)]
    for obj in with_uuid:
        uuid_dict.update({obj.__uuid__: obj})
        uuid_dict.update(get_all_uuids(obj))
    return uuid_dict

def find_dependent_uuids(json_dct):
    dct = ujson.loads(json_dct)
    uuids = [decode_uuid(obj) for obj in flatten_all(dct)
             if is_uuid_string(obj)]
    return uuids

def get_all_uuid_strings(dct):
    all_uuids = []
    for uuid in dct:
        all_uuids.append(uuid)
        all_uuids += find_dependent_uuids(dct[uuid])
    return all_uuids

# NOTE: this only need to find until the first UUID: iterables/mapping with
# UUIDs aren't necessary here
def replace_uuid(obj):
    replacement = obj
    if uuid_test(obj):
        # TODO: compact representation of UUID
        replacement = encode_uuid(obj.__uuid__)
    elif is_mappable(obj):
        replacement = {k: replace_uuid(v) for (k, v) in replacement.items()}
    elif is_iterable(obj):
        replace_type = type(obj)
        replacement = replace_type([replace_uuid(o) for o in obj])
    return replacement

def to_dict_with_uuids(obj):
    dct = obj.to_dict()
    return replace_uuid(dct)

def to_json(obj):
    dct = to_dict_with_uuids(obj)
    dct.update({'__module__': obj.__class__.__module__,
                '__class__': obj.__class__.__name__})
    return ujson.dumps(dct)

def from_json(json_str, existing_uuids):
    # NOTE: from_json only works with existing_uuids
    dct = ujson.loads(json_str)
    mod = importlib.import_module(dct.pop('__module__'))
    cls = getattr(mod, dct.pop('__class__'))
    for (key, value) in dct.items():
        if is_uuid_string(value):
            # raises KeyError if object hasn't been visited
            # (indicates problem in DAG reconstruction)
            value_obj = existing_uuids[decode_uuid(value)]
            dct[key] = value_obj
    return cls.from_dict(dct)

def reconstruction_dag(uuid_json_dict, dag=None)
    dependent_uuids = {uuid: find_dependent_uuids(json_str)
                       for (uuid, json_str) in uuid_json_dict.items()}
    if dag is None:
        dag = nx.DiGraph()
    for from_node, to_nodes in dependent_uuids.items():
        if to_nodes:
            dag.add_edges_from([(from_node, to_node)
                                for to_node in to_nodes])
    if not nx_dag.is_directed_acyclic_graph(dag):
        raise RuntimeError("Reconstruction DAG not acyclic?!?!")
    return dag


# TODO: move this to storage
def serialize(list_of_objects):
    uuid_object_dict = {}
    for obj in list_of_objects:
        uuid_object_dict.update(get_all_uuids(obj))

    # TODO: replace to_json with something that gets the serializer
    uuid_json_dict = {uuid: to_json(obj)
                      for (uuid, obj) in uuid_object_dict.item()}
    return uuid_json_dict


def deserialize(uuid_json_dict, lazies, storage):
    dag = reconstruction_dag(uuid_json_dict)
    missing = check_dag(dag, uuid_json_dict)
    while missing:
        (more_json, loc_lazies) = storage.backend.load_table_data(missing)
        uuid_json_dict.update(more_json)
        lazies = nested_update(lazies, loc_lazies)
        dag = reconstruction_dag(uuid_json_dict, dag)
        missing = check_dag(dag, uuid_json_dict)

    new_uuids = {}
    known_uuids = storage.known_uuids
    for lazy_table in lazies:
        lazy_uuid_objects = {
            lazy.uuid: storage.make_lazy(lazy_table, lazy.uuid)
            for lazy in lazies[lazy_table]
            if lazy.uuid not in known_uuids
        }
        new_uuids.update(lazy_uuid_objects)
        known_uuids.update(lazy_uuid_objects)

    ordered_nodes = list(reversed(list(nx_dag.topological_sort(dag))))
    for node in ordered_nodes:
        # TODO: replace from_json with something that gets the deserializer
        new_uuids[node] = from_json(all_json[node], new_uuids, known_uuids)
    return new_uuids
