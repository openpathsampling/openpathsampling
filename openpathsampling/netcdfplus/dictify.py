import base64
import importlib

import numpy as np

from openpathsampling.integration_tools import is_simtk_quantity, unit
import math
import abc
from uuid import UUID

import ujson

import marshal
import types
import opcode

from .base import StorableObject

from openpathsampling.tools import word_wrap

from .cache import WeakValueCache

__author__ = 'Jan-Hendrik Prinz'

import sys
if sys.version_info > (3, ):
    long = int
    unicode = str
    builtin_module = 'builtins'
    get_code = lambda func: func.__code__
    intify_byte = lambda b: b
    decodebytes = lambda s: base64.decodebytes(s.encode())
    get_callable_name = lambda c: c.__name__
    import builtins
else:
    builtin_module = '__builtin__'
    get_code = lambda func: func.func_code
    intify_byte = lambda b: ord(b)
    decodebytes = base64.decodestring
    get_callable_name = lambda c: c.func_name
    import builtins

# in Python 3.6 the opcodes have changed width
if sys.version_info > (3, 6):
    opcode_arg_width = 1
    opcode_no_arg_width = 1
else:
    opcode_arg_width = 2
    opcode_no_arg_width = 0

if int(ujson.__version__.split(".")[0]) <= 2:
    ujson_kwargs = dict()
else:
    ujson_kwargs = {"reject_bytes": False}

class ObjectJSON(object):
    """
    A simple implementation of a pickle algorithm to create object that can be
    converted to json and back

    Attributes
    ----------
    safemode: bool
        If set to `True` the recreation of marshalled objects like functions is
        switched off and these objects are replaced by None. Can be used to load
        from incompatible python versions or potential unsafe trajectory files.
    """

    allow_marshal = True

    # switch to `true`, if you want more protection
    prevent_unsafe_modules = False

    allowed_storable_atomic_types = [
        int, float, bool, long, str,
        np.float32, np.float64,
        np.int8, np.int16, np.int32, np.int64,
        np.uint8, np.uint16, np.uint32, np.uint64,
    ]

    safe_modules = [
        'numpy',
        'math',
        'pandas',
        'mdtraj',
        'simtk',
        'simtk.unit',
        'simtk.openmm'
    ]

    def __init__(self, unit_system=None):
        self.excluded_keys = []
        self.unit_system = unit_system
        self.class_list = dict()
        self.allowed_storable_types = dict()
        self.type_names = {}
        self.type_classes = {}
        self.safemode = False

        self.update_class_list()

    def update_class_list(self):
        self.class_list = StorableObject.objects()
        self.type_names = {
            cls.__name__: cls for cls in self.allowed_storable_atomic_types}
        self.type_names.update(self.class_list)
        self.type_classes = {
            cls: name for name, cls in self.type_names.items()}

    def simplify_object(self, obj):
        return {
            '_cls': obj.__class__.__name__,
            '_dict': self.simplify(obj.to_dict(), obj.base_cls_name)
        }

    def simplify(self, obj, base_type=''):
        if obj.__class__.__name__ == 'module':
            # store an imported module
            if obj.__name__.split('.')[0] in self.safe_modules:
                return {'_import': obj.__name__}
            else:
                raise RuntimeError((
                    'The module reference "%s" you want to store is '
                    'not allowed!') % obj.__name__)

        elif type(obj) is type or type(obj) is abc.ABCMeta:
            # store a storable number type
            if obj in self.type_classes:
                return {'_type': obj.__name__}
            else:
                return None

        elif type(obj) is float and math.isinf(obj):
            return {
                '_float': str(obj)}

        elif type(obj) is int and math.isinf(obj):
            return {
                '_integer': str(obj)}

        elif obj.__class__.__module__ != builtin_module:
            #if obj.__class__ is units.Quantity:
            if is_simtk_quantity(obj):
                # This is number with a unit so turn it into a list
                if self.unit_system is not None:
                    return {
                        '_value': self.simplify(
                            obj.value_in_unit_system(self.unit_system)),
                        '_units': self.unit_to_dict(
                            obj.unit.in_unit_system(self.unit_system))
                    }
                else:
                    return {
                        '_value': self.simplify(obj / obj.unit, base_type),
                        '_units': self.unit_to_dict(obj.unit)
                    }
            elif obj.__class__ is np.ndarray:
                # this is maybe not the best way to store large numpy arrays!
                return {
                    '_numpy': self.simplify(obj.shape),
                    '_dtype': str(obj.dtype),
                    '_data': base64.b64encode(obj.copy(order='C'))
                }
            elif hasattr(obj, 'to_dict'):
                # the object knows how to dismantle itself into a json string
                if hasattr(obj, '__uuid__'):
                    return {
                        '_cls': obj.__class__.__name__,
                        '_obj_uuid': str(UUID(int=obj.__uuid__)),
                        '_dict': self.simplify(obj.to_dict(), base_type)}
                else:
                    return {
                        '_cls': obj.__class__.__name__,
                        '_dict': self.simplify(obj.to_dict(), base_type)}
            elif type(obj) is UUID:
                return {
                    '_uuid': str(UUID(int=obj))}
            else:
                return None
        elif type(obj) is list:
            return [self.simplify(o, base_type) for o in obj]
        elif type(obj) is tuple:
            return {'_tuple': [self.simplify(o, base_type) for o in obj]}
        elif type(obj) is dict:
            # we want to support storable objects as keys so we need to wrap
            # dicts with care and store them using tuples

            simple = [
                key for key in obj.keys()
                if type(key) is str or type(key) is int]

            if len(simple) < len(obj):
                # other keys than int or str
                result = {
                    '_dict': [
                        self.simplify(tuple([key, o]))
                        for key, o in obj.items()
                        if key not in self.excluded_keys
                    ]}
            else:
                # simple enough, do it the old way
                # FASTER VERSION NORMALLY
                result = {
                    key: self.simplify(o) for key, o in obj.items()
                    if key not in self.excluded_keys
                }

                # SLOWER VERSION FOR DEBUGGING
                # result = {}
                # for key, o in obj.items():
                # logger.debug("Making dict entry of " + str(key) + " : "
                # + str(o))
                # if key not in self.excluded_keys:
                # result[key] = self.simplify(o)
                # else:
                # logger.debug("EXCLUDED")

            return result
        elif type(obj) is slice:
            return {
                '_slice': [obj.start, obj.stop, obj.step]}
        else:
            oo = obj
            return oo

    @staticmethod
    def _unicode2str(s):
        if type(s) is unicode:
            return str(s)
        else:
            return s

    def build(self, obj):
        if type(obj) is dict:
            if '_units' in obj and '_value' in obj:
                return self.build(
                    obj['_value']) * self.unit_from_dict(obj['_units'])

            elif '_slice' in obj:
                return slice(*obj['_slice'])

            elif '_numpy' in obj:
                return np.frombuffer(
                    decodebytes(obj['_data']),
                    dtype=np.dtype(obj['_dtype'])).reshape(
                        self.build(obj['_numpy'])
                )

            elif '_float' in obj:
                return float(str(obj['_float']))

            elif '_integer' in obj:
                return float(str(obj['_integer']))

            elif '_uuid' in obj:
                return int(UUID(obj['_uuid']))

            elif '_cls' in obj and '_dict' in obj:
                if obj['_cls'] not in self.class_list:
                    self.update_class_list()
                    if obj['_cls'] not in self.class_list:
                        # updating did not help, so there is nothing we can do.
                        return None
                        # raise ValueError((
                        #     'Cannot create obj of class `%s`.\n' +
                        #     'Class is not registered as creatable! '
                        #     'You might have to define\n' +
                        #     'the class locally and call '
                        #     '`update_storable_classes()` on your storage.') %
                        #     obj['_cls'])

                attributes = self.build(obj['_dict'])
                return self.class_list[obj['_cls']].from_dict(attributes)

            elif '_tuple' in obj:
                return tuple([self.build(o) for o in obj['_tuple']])

            elif '_type' in obj:
                # return a type of a _built-in_ `netcdfplus` type
                return self.type_names.get(obj['_type'])

            elif '_dict' in obj:
                return {
                    self._unicode2str(self.build(key)): self.build(o)
                    for key, o in self.build(obj['_dict'])
                }

            elif '_import' in obj:
                module = obj['_import']
                if module.split('.')[0] in self.safe_modules:
                    imp = importlib.import_module(module)
                    return imp
                else:
                    return None

            elif '_marshal' in obj or '_module' in obj:
                if self.safemode:
                    return None

                return self.callable_from_dict(obj)

            else:
                return {
                    self._unicode2str(key): self.build(o)
                    for key, o in obj.items()
                }

        elif type(obj) is list:
            return [self.build(o) for o in obj]

        elif type(obj) is unicode:
            return str(obj)

        else:
            return obj

    @staticmethod
    def unit_to_symbol(unit):
        return str(1.0 * unit).split()[1]

    @staticmethod
    def unit_to_dict(unit):
        unit_dict = {
            p.name: int(fac) for p, fac in unit.iter_base_or_scaled_units()}
        return unit_dict

    @staticmethod
    def unit_from_dict(unit_dict):
        # this will *only* work if simtk.unit is installed
        this_unit = unit.Unit({})
        for unit_name, unit_multiplication in unit_dict.items():
            this_unit *= getattr(unit, unit_name) ** unit_multiplication

        return this_unit

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
        root_module = f_module.split('.')[0]

        # is_class = isinstance(c, (type, types.ClassType))

        # try saving known external classes of functions, e.g. `msmbuilder`
        if root_module in ObjectJSON.safe_modules:
            # only store the function/class and the module
            return {
                '_module': c.__module__,
                '_name': c.__name__
            }

        # if the easy way did not work, try saving it using bytecode
        if ObjectJSON.allow_marshal and callable(c):
            # use marshal
            global_vars = ObjectJSON._find_var(c, opcode.opmap['LOAD_GLOBAL'])
            import_vars = ObjectJSON._find_var(c, opcode.opmap['IMPORT_NAME'])

            all_builtins = dir(builtins)

            global_vars = list(set(
                [var for var in global_vars if var not in all_builtins]))
            import_vars = list(set(import_vars))

            err = ''

            if len(global_vars) > 0:
                err += 'The function you try to save relies on globally set ' \
                       'variables and these cannot be saved since storage ' \
                       'has no access to the global scope which includes ' \
                       'imports! \n\n'
                err += 'We require that the following globals: ' + \
                       str(global_vars) + ' either\n'
                err += '\n1. be replaced by constants'
                err += '\n2. be defined inside your function,' + \
                       '\n\n' + '\n'.join(
                           map(lambda x: ' ' * 8 + x + '= ...', global_vars)
                       ) + '\n'
                err += '\n3. imports need to be "re"-imported inside your ' \
                       'function' + \
                       '\n\n' + '\n'.join(
                           map(lambda x: ' ' * 8 + 'import ' + x, global_vars)
                       ) + '\n'
                err += '\n4. be passed as an external parameter ' \
                       '(not for imports!)'
                err += '\n\n        my_cv = FunctionCV("cv_name", ' + \
                       get_callable_name(c) + ', \n' + \
                       ',\n'.join(
                           map(lambda x: ' ' * 20 + x + '=' + x, global_vars)
                       ) + ')' + '\n'
                err += '\n    and change your function definition like this'
                err += '\n\n        def ' + \
                       get_callable_name(c) + '(snapshot, ...,  ' + \
                       '\n' + ',\n'.join(
                           map(lambda x: ' ' * 16 + x, global_vars)
                       ) + '):'

            unsafe_modules = [
                module for module in import_vars
                if module not in ObjectJSON.safe_modules
            ]

            if ObjectJSON.prevent_unsafe_modules and len(unsafe_modules) > 0:
                if len(err) > 0:
                    err += '\n\n'

                err += 'The function you try to save requires the following' \
                       ' modules to be installed: ' + str(unsafe_modules) + \
                       ' which are not marked as safe! '
                err += 'You can change the list of safe modules using '
                err += '\n\n        ObjectJSON.safe_modules.extend(['
                err += '\n' + ',\n'.join(
                       map(lambda x: ' ' * 12 + x, unsafe_modules)
                )
                err += '\n        ])'
                err += '\n\n'
                err += 'include the import statement in your function like'
                err += '\n\n' + '\n'.join(
                    [' ' * 8 + 'import ' + v for v in unsafe_modules])

            if len(err) > 0:
                raise RuntimeError('Cannot store function! \n\n' +
                                   word_wrap(err, 60))

            return {
                '_marshal': ObjectJSON._to_marshal(c),
                '_global_vars': global_vars,
                '_module_vars': import_vars
            }

        raise RuntimeError('Locally defined classes are not storable yet')

    @staticmethod
    def _to_marshal(c):
        return base64.b64encode(marshal.dumps(get_code(c)))

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
                if packages[0] in ObjectJSON.safe_modules:
                    imp = importlib.import_module(module)
                    c = getattr(imp, c_dict['_name'])

        return c

    @staticmethod
    def _to_opcode(code):
        """
       Yields the tuple opcode + argument for code

        Parameters
        ----------
        code : function
            the python bytecode to be searched

        Returns
        -------
        generator of (int, int)
            generator the returns tuples of opcode + argument
        """

        opcode_stream = get_code(code).co_code
        i = 0
        extended = 0
        length = len(opcode_stream)
        while i < length:
            int_code = intify_byte(opcode_stream[i])
            i += 1
            if int_code >= opcode.HAVE_ARGUMENT:
                opargs = intify_byte(opcode_stream[i]) + extended
                if opcode_arg_width == 2:
                    # before Python 3.6
                    opargs += intify_byte(opcode_stream[i + 1]) * 256

                extended = 0
                i += opcode_arg_width

                if int_code == opcode.EXTENDED_ARG:
                    extended = (256 ** opcode_arg_width) * opargs
                    continue

            else:
                i += opcode_no_arg_width
                opargs = None

            yield int_code, opargs

    @staticmethod
    def _find_var(code, op):
        """
        Helper function to search in python bytecode for specific function calls

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

        code_object = get_code(code)
        if op in opcode.hasconst:
            variable = code_object.co_consts
        elif op in opcode.haslocal:
            variable = code_object.co_varnames
        elif op in opcode.hasname:
            variable = code_object.co_names
        elif op in opcode.hasfree:
            variable = code_object.co_freevars
        else:
            return []

        return list(set([
            variable[arg] for code, arg in ObjectJSON._to_opcode(code)
            if code == op and arg is not None]))

    def to_json(self, obj, base_type=''):
        simplified = self.simplify(obj, base_type)
        return ujson.dumps(simplified, **ujson_kwargs)

    def to_json_object(self, obj):
        if hasattr(obj, 'base_cls') \
                and type(obj) is not type and type(obj) is not abc.ABCMeta:
            simplified = self.simplify_object(obj)
        else:
            simplified = self.simplify(obj)
        try:
            json_str = ujson.dumps(simplified, **ujson_kwargs)
        except TypeError as e:
            err = (
                'Cannot convert object of type `%s` to json. '
                '\n__dict__: %s\n'
                '\nsimplified: %s\n'
                '\nError: %s'
            ) % (
                obj.__class__.__name__,
                obj.__dict__,
                simplified,
                str(e)
            )
            raise ValueError(err)

        return json_str

    def from_json(self, json_string):
        simplified = ujson.loads(json_string)
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
        if obj.__class__.__module__ != builtin_module:
            if obj.__class__ in self.storage._obj_store:
                store = self.storage._obj_store[obj.__class__]
                if not store.nestable or obj.base_cls_name != base_type:
                    # this also returns the base class name used for storage
                    # store objects only if they are not creatable. If so they
                    # will only be created in their top instance and we use
                    # the simplify from the super class ObjectJSON
                    idx = store.save(obj)
                    if idx is None:
                        raise RuntimeError(
                            'cannot store idx None in store %s' % store)
                    return {
                        '_idx': idx,
                        '_store': store.prefix}

        return super(StorableObjectJSON, self).simplify(obj, base_type)

    def build(self, obj):
        if type(obj) is dict:
            if '_storage' in obj:
                if obj['_storage'] == 'self':
                    return self.storage

            if '_idx' in obj and '_store' in obj:
                store = self.storage._stores[obj['_store']]
                result = store.load(obj['_idx'])

                return result

        return super(StorableObjectJSON, self).build(obj)


class UUIDObjectJSON(ObjectJSON):
    def __init__(self, storage, unit_system=None):
        super(UUIDObjectJSON, self).__init__(unit_system)
        self.excluded_keys = ['json']
        self.storage = storage

    def simplify(self, obj, base_type=''):
        if obj is self.storage:
            return {'_storage': 'self'}

        if obj.__class__.__module__ != builtin_module:
            if obj.__class__ in self.storage._obj_store:
                store = self.storage._obj_store[obj.__class__]
                if not store.nestable or obj.base_cls_name != base_type:
                    # this also returns the base class name used for storage
                    # store objects only if they are not creatable. If so
                    # they will only be created in their top instance and we
                    # use the simplify from the super class ObjectJSON
                    store.save(obj)
                    return {
                        '_hex_uuid': hex(obj.__uuid__),
                        '_store': store.prefix}
                    # return {
                    #     '_obj_uuid': str(UUID(int=obj.__uuid__)),
                    #     '_store': store.prefix}

        return super(UUIDObjectJSON, self).simplify(obj, base_type)

    def build(self, obj):
        if type(obj) is dict:
            if '_storage' in obj:
                if obj['_storage'] == 'self':
                    return self.storage

            if '_obj_uuid' in obj and '_store' in obj:
                store = self.storage._stores[obj['_store']]
                result = store.load(int(UUID(obj['_obj_uuid'])))

                return result

            if '_hex_uuid' in obj and '_store' in obj:
                store = self.storage._stores[obj['_store']]
                result = store.load(int(obj['_hex_uuid'].strip('L'), 16))

                return result

        return super(UUIDObjectJSON, self).build(obj)


class CachedUUIDObjectJSON(ObjectJSON):
    def __init__(self, unit_system=None):
        super(CachedUUIDObjectJSON, self).__init__(unit_system)
        self.excluded_keys = ['json']
        self.uuid_cache = WeakValueCache()

    def simplify(self, obj, base_type=''):
        if obj.__class__.__module__ != builtin_module:
            if hasattr(obj, 'to_dict') and hasattr(obj, '__uuid__'):
                # the object knows how to dismantle itself into a json string
                if obj.__uuid__ not in self.uuid_cache:
                    self.uuid_cache[obj.__uuid__] = obj

                    return {
                        '_cls': obj.__class__.__name__,
                        '_obj_uuid': str(UUID(int=obj.__uuid__)),
                        '_dict': self.simplify(obj.to_dict(), base_type)}
                else:
                    return {
                        '_obj_uuid': str(UUID(int=obj.__uuid__))}

        return super(CachedUUIDObjectJSON, self).simplify(obj, base_type)

    def build(self, jsn):
        if type(jsn) is dict:
            if '_obj_uuid' in jsn:
                uuid = UUID(jsn['_obj_uuid'])
                if uuid in self.uuid_cache:
                    return self.uuid_cache[uuid]
                elif '_cls' in jsn and '_dict' in jsn:
                    if jsn['_cls'] not in self.class_list:
                        self.update_class_list()
                        if jsn['_cls'] not in self.class_list:
                            raise ValueError((
                                 'Cannot create jsn of class `%s`.\n' +
                                 'Class is not registered as creatable! '
                                 'You might have to define\n' +
                                 'the class locally and call '
                                 '`update_storable_classes()` on your storage.'
                            ) % jsn['_cls'])

                    attributes = self.build(jsn['_dict'])

                    obj = self.class_list[jsn['_cls']].from_dict(attributes)
                    obj.__uuid__ = uuid
                    self.uuid_cache[uuid] = obj
                    return obj
                else:
                    # this should not happen!
                    raise ('What happend here. JSN `%s`' % jsn)
                    pass

        return super(CachedUUIDObjectJSON, self).build(jsn)

    def to_json(self, obj, base_type=''):
        # we need to clear the cache, since we have no idea, what the other end
        # still knows. We can only cache stuff we are sending this time
        self.uuid_cache.clear()
        return super(CachedUUIDObjectJSON, self).to_json(obj, base_type)

    # def from_json(self, json_string):
    #     # here we keep the cache. It could happen that an object is sent in
    #     # full, but we still have it and so we do not have to rebuild it which
    #     # saves some time
    #     simplified = ujson.loads(json_string)
    #     return self.build(simplified)
