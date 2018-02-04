import collections

import sys
if sys.version_info > (3, ):
    basestring = str
else:
    basestring = basestring

def is_mappable(obj):
    return isinstance(obj, collections.Mapping)

def is_iterable(obj):
    return (isinstance(obj, collections.Iterable)
            and not isinstance(obj, basestring))


def none_to_default(option, default):
    if option is None:
        option = default
    return option

def group_by(list_of_iterable, group, grouping_function):
    results = collections.defaultdict(list)
    for obj in list_of_iterable:
        results[grouping_function(obj, group)].append(obj)
    return results

def group_by_index(list_of_iterable, column_number):
    results = collections.defaultdict(list)
    for obj in list_of_iterable:
        results[obj[column_number]].append(obj)
    return results

def group_by_attribute(list_of_iterable, attr):
    return group_by(list_of_iterable, attr, getattr)

def dict_group_by(dct, group, key_extract):
    results = collections.defaultdict(dict)
    for (key, value) in dct.items():
        results[key_extract(key, value)].update({key: value})
    return results


def compare_sets(set1, set2):
    only_in_1 = set1 - set2
    only_in_2 = set2 - set1
    return (only_in_1, only_in_2)

def flatten(inputs, value_iter, classes):
    results = []
    for val in value_iter(inputs):
        if isinstance(val, classes) and not isinstance(val, basestring):
            results += flatten(val, value_iter, classes)
        else:
            results.append(val)
    return results

def flatten_dict(dct):
    return flatten(dct, lambda x: x.values(), dict)

def flatten_iterable(ll):
    return flatten(ll, lambda x: x, (list, tuple, set))

def flatten_all(obj):
    is_mappable = lambda x: isinstance(x, collections.Mapping)
    return flatten(obj,
                   lambda x: x.values() if is_mappable(x) else x.__iter__(),
                   (collections.Mapping, collections.Iterable))

def nested_update(original, update):
    for (k, v) in update.items():
        if isinstance(v, collections.Mapping):
            original[k] = nested_update(original.get(k, {}), v)
        else:
            original[k] = v
    return original

