import itertools
import collections
from collections import abc
from numpy import ndarray

import logging
logger = logging.getLogger(__name__)

class SimpleNamespace(abc.MutableMapping):
    # types.SimpleNameSpace in 3.3+
    # this variants acts as either a dict or a namespace for getting
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __delitem__(self, item):
        del self.__dict__[item]

    def __len__(self):
        return len(self.__dict__)

    def __iter__(self):
        return self.__dict__.__iter__()

    def __repr__(self):
        keys = sorted(self.__dict__)
        items = ("{}={!r}".format(k, self.__dict__[k]) for k in keys)
        return "{}({})".format(type(self).__name__, ", ".join(items))

    def __eq__(self, other):
        return self.__dict__ == other.__dict__


import sys
if sys.version_info > (3, ):
    basestring = str
else:
    basestring = basestring

# simplifications for the necessary type-checking
def is_string(obj):
    return isinstance(obj, basestring)

def is_mappable(obj):
    return isinstance(obj, abc.Mapping)

def is_iterable(obj):
    return isinstance(obj, abc.Iterable) and not is_string(obj)

def is_numpy_iterable(obj):
    return isinstance(obj, ndarray)

def listify(obj):
    if not is_iterable(obj):
        obj = [obj]
    return obj

def none_to_default(option, default):
    if option is None:
        option = default
    return option

# group_by and variants
def group_by(list_of_iterable, group, grouping_function):
    results = collections.defaultdict(list)
    for obj in list_of_iterable:
        results[grouping_function(obj, group)].append(obj)
    return results

def group_by_function(ll, function):
    results = collections.defaultdict(list)
    for obj in ll:
        results[function(obj)].append(obj)
    return results

def group_by_index(list_of_iterable, column_number):
    results = collections.defaultdict(list)
    for obj in list_of_iterable:
        results[obj[column_number]].append(obj)
    return results

def group_by_attribute(list_of_iterable, attr):
    return group_by(list_of_iterable, attr, getattr)

def dict_group_by(dct, key_extract):
    results = collections.defaultdict(dict)
    for (key, value) in dct.items():
        results[key_extract(key, value)].update({key: value})
    return results

def block(sliceable, length):
    sliceable = list(sliceable)
    max_len = len(sliceable)
    n_iterations = max_len // length + 1
    min_val = 0
    max_val = min(length, max_len)
    while min_val <= max_len:
        yield sliceable[slice(min_val, max_val)]
        min_val += length
        max_val += length
        max_val = min(max_val, max_len)

_no_fill = object()
def grouper(iterable, length, fillvalue=_no_fill):
    # based on https://stackoverflow.com/a/10791887
    args = [iter(iterable)] * length
    for block in itertools.zip_longest(*args, fillvalue=fillvalue):
        if tuple(block)[-1] is _no_fill:
            yield tuple(elem for elem in block if elem is not _no_fill)
        else:
            yield tuple(block)

def compare_sets(set1, set2):
    only_in_1 = set1 - set2
    only_in_2 = set2 - set1
    return (only_in_1, only_in_2)

# flatten and variants
def flatten(inputs, value_iter, classes, excluded=None):
    excluded = none_to_default(excluded, basestring)
    results = []
    for val in value_iter(inputs):
        if isinstance(val, classes) and not isinstance(val, excluded):
            results += flatten(val, value_iter, classes, excluded)
        else:
            results.append(val)
    return results


def flatten_dict(dct):
    return flatten(dct, lambda x: x.values(), dict)

def flatten_iterable(ll):
    return flatten(ll, lambda x: x, (list, tuple, set))

def flatten_all(obj):
    is_mappable = lambda x: isinstance(x, abc.Mapping)
    return flatten(obj,
                   lambda x: x.values() if is_mappable(x) else x.__iter__(),
                   (abc.Mapping, abc.Iterable),
                   (basestring, ndarray))

def nested_update(original, update):
    for (k, v) in update.items():
        if isinstance(v, abc.Mapping):
            original[k] = nested_update(original.get(k, {}), v)
        else:
            original[k] = v
    return original

