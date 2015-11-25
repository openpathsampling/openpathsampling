import base64
import json
import importlib

import numpy as np
from simtk import unit as units
import yaml
import abc

from base import StorableObject

__author__ = 'Jan-Hendrik Prinz'


class ObjectJSON(object):
    """
    A simple implementation of a pickle algorithm to create object that can be converted to json and back
    """

    allow_marshal = True

    allowed_storable_atomic_types = [
        int, float, bool, long, str,
        np.float32, np.float64,
        np.int8, np.int16, np.int32, np.int64,
        np.uint8, np.uint16, np.uint32, np.uint64,
    ]

    allowed_imports = [
        'numpy',
        'math',
        'pandas',
        'mdtraj'
    ]

    def __init__(self, unit_system=None):
        self.excluded_keys = []
        self.unit_system = unit_system
        self.class_list = dict()
        self.allowed_storable_types = dict()

        self.update_class_list()

    def update_class_list(self):
        self.class_list = StorableObject.objects()
        self.type_names = {cls.__name__: cls for cls in self.allowed_storable_atomic_types}
        self.type_names.update(self.class_list)
        self.type_classes = {cls: name for name, cls in self.type_names.iteritems()}

    def simplify_object(self, obj):
        return {'_cls': obj.__class__.__name__, '_dict': self.simplify(obj.to_dict(), obj.base_cls_name)}

    def simplify(self, obj, base_type=''):
        if obj.__class__.__name__ == 'module':
            # store an imported module
            if obj.__name__.split('.')[0] in self.allowed_imports:
                return {'_import': obj.__name__}
            else:
                raise RuntimeError('The module reference "%s" you want to store is not allowed!' % obj.__name__)
        elif type(obj) is type or type(obj) is abc.ABCMeta:
            # store a storable number type
            if obj in self.type_classes:
                return {'_type': obj.__name__}
            else:
                return None

        elif obj.__class__.__module__ != '__builtin__':
            if obj.__class__ is units.Quantity:
                # This is number with a unit so turn it into a list
                if self.unit_system is not None:
                    return {'_value': obj.value_in_unit_system(self.unit_system),
                            '_units': self.unit_to_dict(obj.unit.in_unit_system(self.unit_system))}
                else:
                    return {'_value': obj / obj.unit, '_units': self.unit_to_dict(obj.unit)}
            elif obj.__class__ is np.ndarray:
                # this is maybe not the best way to store large numpy arrays!
                return {'_numpy': self.simplify(obj.shape), '_dtype': str(obj.dtype), '_data': base64.b64encode(obj)}
            elif hasattr(obj, 'to_dict'):
                # the object knows how to dismantle itself into a json string so use this
                return {'_cls': obj.__class__.__name__, '_dict': self.simplify(obj.to_dict(), base_type)}
            else:
                return None
        elif type(obj) is list:
            return [self.simplify(o, base_type) for o in obj]
        elif type(obj) is tuple:
            return {'_tuple': [self.simplify(o, base_type) for o in obj]}
        elif type(obj) is dict:
            # we want to support storable objects as keys so we need to wrap
            # dicts with care and store them using tuples

            simple = [key for key in obj.keys() if type(key) is str or type(key) is int]

            if len(simple) < len(obj):
                # other keys than int or str
                result = {'_dict':
                              [
                                self.simplify(tuple([key, o]))
                                for key, o in obj.iteritems()
                                if key not in self.excluded_keys
                              ]
                          }

            else:
                # simple enough, do it the old way
                # FASTER VERSION NORMALLY
                result = {key: self.simplify(o)
                          for key, o in obj.iteritems()
                          if key not in self.excluded_keys
                          }

                # SLOWER VERSION FOR DEBUGGING
                # result = {}
                # for key, o in obj.iteritems():
                # logger.debug("Making dict entry of " + str(key) + " : "
                # + str(o))
                # if key not in self.excluded_keys:
                # result[key] = self.simplify(o)
                # else:
                # logger.debug("EXCLUDED")

            return result
        elif type(obj) is slice:
            return {'_slice': [obj.start, obj.stop, obj.step]}
        else:
            oo = obj
            return oo

    def build(self, obj):
        if type(obj) is dict:
            if '_units' in obj and '_value' in obj:
                return obj['_value'] * self.unit_from_dict(obj['_units'])
            elif '_slice' in obj:
                return slice(*obj['_slice'])
            elif '_numpy' in obj:
                return np.frombuffer(base64.decodestring(obj['_data']), dtype=np.dtype(obj['_dtype'])).reshape(
                    self.build(obj['_numpy']))
            elif '_cls' in obj and '_dict' in obj:
                if obj['_cls'] not in self.class_list:
                    self.update_class_list()
                    if obj['_cls'] not in self.class_list:
                        # updating did not help, so there is nothing we can do.
                        raise ValueError('Cannot create obj of class "' + obj['_cls'] + '".\n' +
                                         'Class is not registered as creatable! You might have to define\n' +
                                         'the class locally and call update_storable_classes() on your storage.')

                attributes = self.build(obj['_dict'])
                return self.class_list[obj['_cls']].from_dict(attributes)
            elif '_tuple' in obj:
                return tuple([self.build(o) for o in obj['_tuple']])
            elif '_type' in obj:
                # return a type of a built-in type that represents a type in netcdf
                return self.type_names.get(obj['_type'])

            elif '_dict' in obj:
                return {
                    self.build(key): self.build(o)
                    for key, o in self.build(obj['_dict'])
                    }
            elif '_import' in obj:
                module = obj['_import']
                if module.split('.')[0] in self.allowed_imports:
                    imp = importlib.import_module(module)
                    return imp
                else:
                    return None

            else:
                return {
                    key: self.build(o)
                    for key, o in obj.iteritems()
                    }
        elif type(obj) is list:
            return [self.build(o) for o in obj]
        else:
            return obj

    @staticmethod
    def unitsytem_to_list(unit_system):
        """
        Turn a simtk.UnitSystem() into a list of strings representing the unitsystem for serialization
        """
        return [u.name for u in unit_system.units]

    @staticmethod
    def unit_system_from_list(unit_system_list):
        """
        Create a simtk.UnitSystem() from a serialialized list of strings representing the unitsystem
        """
        return units.UnitSystem(
            [getattr(units, unit_name).iter_base_or_scaled_units().next()[0] for unit_name in unit_system_list])

    @staticmethod
    def unit_to_symbol(unit):
        return str(1.0 * unit).split()[1]

    @staticmethod
    def unit_to_dict(unit):
        unit_dict = {p.name: int(fac) for p, fac in unit.iter_base_or_scaled_units()}
        return unit_dict

    @staticmethod
    def unit_from_dict(unit_dict):
        unit = units.Unit({})
        for unit_name, unit_multiplication in unit_dict.iteritems():
            unit *= getattr(units, unit_name) ** unit_multiplication

        return unit

    def to_json(self, obj, base_type=''):
        simplified = self.simplify(obj, base_type)
        return json.dumps(simplified)

    def to_json_object(self, obj):
        if hasattr(obj, 'base_cls') and type(obj) is not type and type(obj) is not abc.ABCMeta:
            simplified = self.simplify_object(obj)
        else:
            simplified = self.simplify(obj)
        try:
            json_str = json.dumps(simplified)
        except TypeError:
            print obj.__class__.__name__
            print obj.__dict__
            print simplified
            raise ValueError('Not possible to turn object into json')

        return json_str

    def from_json(self, json_string):
        simplified = yaml.load(json_string)
        return self.build(simplified)

    def unit_to_json(self, unit):
        simple = self.unit_to_dict(unit)
        return self.to_json(simple)

    def unit_from_json(self, json_string):
        return self.unit_from_dict(self.from_json(json_string))


class StorableObjectJSON(ObjectJSON):
    def __init__(self, storage, unit_system=None):
        super(StorableObjectJSON, self).__init__(unit_system)
        self.excluded_keys = ['idx', 'json', 'identifier']
        self.storage = storage

    def simplify(self, obj, base_type=''):
        if obj is self.storage:
            return {'_storage': 'self'}
        if obj.__class__.__module__ != '__builtin__':
            if obj.__class__ in self.storage._obj_store:
                store = self.storage._obj_store[obj.__class__]
                if not store.nestable or obj.base_cls_name != base_type:
                    # this also returns the base class name used for storage
                    # store objects only if they are not creatable. If so they will only be created in their
                    # top instance and we use the simplify from the super class ObjectJSON
                    self.storage.save(obj)
                    return {'_idx': store.index[obj], '_obj': store.prefix}

        return super(StorableObjectJSON, self).simplify(obj, base_type)

    def build(self, obj):
        if type(obj) is dict:
            if '_storage' in obj:
                if obj['_storage'] == 'self':
                    return self.storage

            if '_idx' in obj and '_obj' in obj:
                store = self.storage._objects[obj['_obj']]
                result = store.load(obj['_idx'])

                return result

        return super(StorableObjectJSON, self).build(obj)
