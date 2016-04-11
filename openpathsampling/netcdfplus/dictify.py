import base64
import json
import importlib

import numpy as np
from simtk import unit as units
import yaml
import abc

import marshal
import types
import opcode
import __builtin__

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

            elif '_marshal' in obj or '_module' in obj:
                return self.callable_from_dict(obj)

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

    @staticmethod
    def callable_to_dict(c):
        """
        Turn a callable function of class into a dictionary

        Used for conversion to JSON

        Parameters
        ----------
        c : callable (function or class with __call__)
            the function to be turned into a dict representation

        Returns
        -------
        dict
            the dict representation of the callable
        """
        f_module = c.__module__
        is_local = f_module == '__main__'
        is_loaded = f_module.split('.')[0] == 'openpathsampling'
        is_class = isinstance(c, (type, types.ClassType))
        if not is_class:
            if is_local or is_loaded:
                # this is a local function, let's see if we can save it
                if ObjectJSON.allow_marshal and callable(c):
                    # use marshal
                    global_vars = ObjectJSON._find_var(c, opcode.opmap['LOAD_GLOBAL'])
                    import_vars = ObjectJSON._find_var(c, opcode.opmap['IMPORT_NAME'])

                    builtins = dir(__builtin__)

                    global_vars = list(set([
                                               var for var in global_vars if var not in builtins
                                               ]))

                    import_vars = list(set(import_vars))

                    if len(global_vars) > 0:
                        err = 'The function you try to save relies on globally set variables' + \
                              '\nand these cannot be saved since storage has no access to the' + \
                              '\nglobal scope. This includes imports!'
                        err += '\nWe require that the following globals: ' + str(global_vars) + ' either'
                        err += '\n(1) be replaced by constants'
                        err += '\n(2) be defined inside your function,' + \
                               '\n' + '\n'.join(map(lambda x: '    ' + x + '= ...', global_vars))
                        err += '\n(3) imports need to be "re"-imported inside your function' + \
                               '\n' + '\n'.join(map(lambda x: '    import ' + x, global_vars))
                        err += '\n(4) be passed as an external parameter (does not work for imports!), like in '
                        err += '\n        my_cv = CV_Function("cv_name", ' + c.func_name + ', ' + \
                               ', '.join(map(lambda x: x + '=' + x, global_vars)) + ')'
                        err += '\n    and change your function definition like this'
                        err += '\n        def ' + c.func_name + '(snapshot, ...,  ' + \
                               ', '.join(global_vars) + '):'

                        print err

                        raise RuntimeError('Cannot store function! Dependency on global variables')

                        # print [obj._idx for obj in global_vars if hasattr(obj, '_idx')]
                        # print [obj for obj in global_vars]

                    not_allowed_modules = [module for module in import_vars
                                           if module not in ObjectJSON.allowed_imports]

                    if len(not_allowed_modules) > 0:
                        err = 'The function you try to save requires the following modules to ' + \
                              '\nbe installed: ' + str(not_allowed_modules) + ' which are not marked as safe!'
                        err += '\nYou can change the list of safe modules in "CV_function._allowed_modules"'
                        err += '\nYou can also include the import startement in your function like'
                        err += '\n' + '\n'.join(['import ' + v for v in not_allowed_modules])

                        print err

                        raise RuntimeError('Cannot store function! Not allowed modules used.')

                    return {
                        '_marshal': base64.b64encode(
                            marshal.dumps(c.func_code)),
                        '_global_vars': global_vars,
                        '_module_vars': import_vars
                    }

        if not is_local:
            # save the external class, e.g. msmbuilder featurizer
            if f_module.split('.')[0] in ObjectJSON.allowed_imports:
                # only store the function and the module
                return {
                    '_module': c.__module__,
                    '_name': c.__name__
                }

        raise RuntimeError('Locally defined classes are not storable yet')

    @staticmethod
    def callable_from_dict(c_dict):
        """
        Turn a dictionary back in a callable function or class

        Used for conversion from JSON

        Parameters
        ----------
        c_dict : dict
            the dictionary that contains the information

        Returns
        -------
        callable
            the reconstructed callable function or class

        """
        c = None

        if c_dict is not None:
            if '_marshal' in c_dict:
                if ObjectJSON.allow_marshal:
                    code = marshal.loads(base64.b64decode(c_dict['_marshal']))
                    c = types.FunctionType(code, globals(), code.co_name)

            elif '_module' in c_dict:
                module = c_dict['_module']
                packages = module.split('.')
                if packages[0] in ObjectJSON.allowed_imports:
                    imp = importlib.import_module(module)
                    c = getattr(imp, c_dict['_name'])

        return c

    @staticmethod
    def _find_var(code, op):
        """
        Helper function to search in python bytecode for a specific function call

        Parameters
        ----------
        code : function
            the python bytecode to be searched
        op : int
            the int code of the code to be found

        Returns
        -------
        list of func_code.co_names
            a list of co_names used in this function when calling op
        """

        #TODO: Clean this up. It now works only for codes that use co_names
        opcodes = code.func_code.co_code
        i = 0
        ret = []
        while i < len(opcodes):
            int_code = ord(opcodes[i])
            if int_code == op:
                ret.append((i, ord(opcodes[i + 1]) + ord(opcodes[i + 2]) * 256))

            if int_code < opcode.HAVE_ARGUMENT:
                i += 1
            else:
                i += 3

        return [code.func_code.co_names[i[1]] for i in ret]

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
                store = self.storage._stores[obj['_obj']]
                result = store.load(obj['_idx'])

                return result

        return super(StorableObjectJSON, self).build(obj)
