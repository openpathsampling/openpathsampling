import base64
import json
import mdtraj as md
import numpy as np
import pandas as pd
from simtk import unit as units
import simtk.openmm
import yaml
import openpathsampling as paths

class ObjectJSON(object):
    """
    A simple implementation of a pickle algorithm to create object that can be converted to json and back
    """

    def __init__(self, unit_system = None, class_list = None):
        self.excluded_keys = ['name']
        self.unit_system = unit_system
        if class_list is not None:
            self.class_list = class_list
        else:
            self.class_list = paths.todict.class_list

    def simplify_object(self, obj, base_type = ''):
        return { '_cls' : obj.__class__.__name__, '_dict' : self.simplify(obj.to_dict(), obj.base_cls_name) }

    def simplify(self, obj, base_type = ''):
        if type(obj).__module__ != '__builtin__':
            if type(obj) is units.Quantity:
                # This is number with a unit so turn it into a list
                if self.unit_system is not None:
                    return { '_value' : obj.value_in_unit_system(self.unit_system), '_units' : self.unit_to_dict(obj.unit.in_unit_system(self.unit_system)) }
                else:
                    return { '_value' : obj / obj.unit, '_units' : self.unit_to_dict(obj.unit) }
            elif type(obj) is np.ndarray:
                # this is maybe not the best way to store large numpy arrays!
                return { '_numpy' : self.simplify(obj.shape), '_dtype' : str(obj.dtype), '_data' : base64.b64encode(obj) }
            elif hasattr(obj, 'to_dict'):
                # the object knows how to dismantle itself into a json string so use this
                return { '_cls' : obj.__class__.__name__, '_dict' : self.simplify(obj.to_dict(), base_type) }
            else:
                return None
        elif type(obj) is list:
            return [self.simplify(o, base_type) for o in obj]
        elif type(obj) is tuple:
            return tuple([self.simplify(o, base_type) for o in obj])
        elif type(obj) is dict:
            result = {key : self.simplify(o, base_type) for key, o in obj.iteritems() if type(key) is str and key not in self.excluded_keys }
            return result
        elif type(obj) is slice:
            return { '_slice' : [obj.start, obj.stop, obj.step]}
        else:
            oo = obj
            return oo

    def build(self,obj):
        global class_list
        if type(obj) is dict:
            if '_units' in obj and '_value' in obj:
                return obj['_value'] * self.unit_from_dict(obj['_units'])
            elif '_slice' in obj:
                return slice(*obj['_slice'])
            elif '_numpy' in obj:
                return np.frombuffer(base64.decodestring(obj['_data']), dtype=np.dtype(obj['_dtype'])).reshape(tuple(obj['_numpy']))
            elif '_cls' in obj and '_dict' in obj:
                if obj['_cls'] in self.class_list:
                    attributes = self.build(obj['_dict'])
                    return self.class_list[obj['_cls']].from_dict(attributes)
                else:
                    raise ValueError('Cannot create obj of class "' + obj['_cls']+ '". Class is not registered as creatable!')
            else:
                return {key : self.build(o) for key, o in obj.iteritems()}
        elif type(obj) is tuple:
            return tuple([self.build(o) for o in obj])
        elif type(obj) is list:
            return [self.build(o) for o in obj]
        else:
            return obj

    def unitsytem_to_list(self, unit_system):
        '''
        Turn a simtk.UnitSystem() into a list of strings representing the unitsystem for serialization
        '''
        return [ u.name  for u in unit_system.units ]

    def unit_system_from_list(self, unit_system_list):
        '''
        Create a simtk.UnitSystem() from a serialialized list of strings representing the unitsystem
        '''
        return units.UnitSystem([ getattr(units, unit_name).iter_base_or_scaled_units().next()[0] for unit_name in unit_system_list])

    def unit_to_symbol(self, unit):
        return str(1.0 * unit).split()[1]

    def unit_to_dict(self, unit):
        unit_dict = {p.name : int(fac) for p, fac in unit.iter_base_or_scaled_units()}
        return unit_dict

    def unit_from_dict(self, unit_dict):
        unit = units.Unit({})
        for unit_name, unit_multiplication in unit_dict.iteritems():
            unit *= getattr(units, unit_name)**unit_multiplication

        return unit

    def to_json(self, obj, base_type = ''):
        simplified = self.simplify(obj, base_type)
        return json.dumps(simplified)

    def to_json_object(self, obj, base_type = ''):
        simplified = self.simplify_object(obj, base_type)
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

    def topology_to_json(self, topology):
        return self.to_json(self.topology_to_dict(topology))

    def topology_from_json(self, json_string):
        return self.topology_from_dict(self.from_json(json_string))

    def unit_to_json(self, unit):
        simple = self.unit_to_dict(unit)
        return self.to_json(simple)

    def unit_from_json(self, json_string):
        return self.unit_from_dict(self.from_json(json_string))

# Register a class to be creatable. Basically just a dict to match a classname to the actual class
# This is mainly for security so that we do not have to use globals to find classes

class_list = dict()

def restores_as_full_object(super_class):
    class_list[super_class.__name__] = super_class

    super_class.creatable = True
    if not hasattr(super_class, '_excluded_attr'):
        super_class._excluded_attr = []

    if not hasattr(super_class, 'to_dict'):
        def _to_dict(self):
            excluded_keys = ['idx']
            return {key: value for key, value in self.__dict__.iteritems() if key not in excluded_keys and key not in self._excluded_attr and not key.startswith('_')}

        super_class.to_dict = _to_dict

    if not hasattr(super_class, 'from_dict'):
        def _from_dict(cls, my_dict = None):
            if my_dict is None:
                my_dict={}
            try:
                obj = cls(**my_dict)
            except TypeError as e:
                print my_dict
                print cls.__name__
                print e
            return obj

        super_class.from_dict = classmethod(_from_dict)

    return super_class


class LoadedObject(object):
    @property
    def cls(self):
        return self._cls


def restores_as_stub_object(super_class):
    """
    A class decorator that marks a class to be storable in the storage using a LoadedObject class.
    This object will have the same class name, the same dict, but none of the functions and will not have
    been initialized. If you want real objects use @creatable
    :param super_class: The class to be decorated
    :return: The decorated class
    """
    class_list[super_class.__name__] = super_class

    super_class.dictable = True
    if not hasattr(super_class, '_excluded_attr'):
        super_class._excluded_attr = []

    if not hasattr(super_class, 'to_dict'):
        def _to_dict(self):
            excluded_keys = ['idx']
            return {key: value for key, value in self.__dict__.iteritems() if key not in excluded_keys and key not in self._excluded_attr and not key.startswith('_')}

        super_class.to_dict = _to_dict

    if not hasattr(super_class, 'from_dict'):
        def _from_dict(cls, my_dict = None):
            if my_dict is None:
                my_dict={}

            obj = LoadedObject()

            for key, value in my_dict.iteritems():
                setattr(obj, key, value)

            setattr(obj, '_cls', cls.__name__)

            return obj

        super_class.from_dict = classmethod(_from_dict)

    return super_class