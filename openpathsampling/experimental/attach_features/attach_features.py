from collections import namedtuple
import functools
from copy import deepcopy
import inspect
from numbers import Number
import numpy as np

from openpathsampling.netcdfplus import StorableObject

from .write_code import make_init, make_copy_with_replacement

def set_uuid(obj, uuid):
    obj.__uuid__ = uuid

def get_uuid(obj):
    return obj.__uuid__

# we use inpect.Parameter.empty in order to more easily convert to
# inspect.Parameter
_empty = inspect.Parameter.empty

class Parameter:
    """Snapshot parameter description.

    Parameters
    ----------
    name: str
        Parameter name
    default: Union[Int, Float Str, None]
        Default value for this parameter. Must be numeric, string, or None.
        If no default value, use the special value Parameter.empty.
    lazy: bool
        Whether this parameter should use lazy loading
    copy: Callable[[Snapshot], Snapshot]
        Method used to copy this function. Default ``None`` gives
        copy.deepcopy.
    operations: Dict[str, Callable[[Any], None]
        Custom operations that can be performed on this parameter. Keys are
        the method name, values are callables that change this parameter
        based
    docstring: str
        The docstring to show for this parmaeter. This should be in numpydoc
        format, and should include the paramter name on the first line.
    """
    empty = _empty
    def __init__(self, name, default=_empty, lazy=False, copy=None,
                 operations=None, docstring=None):
        self.name = name
        self.default = default
        self.lazy = lazy
        if copy is None:
            copy = deepcopy
        self.copy = copy
        if operations is None:
            operations = {}
        self.operations = operations
        for op_name, op in operations.items():
            setattr(self, op_name, op)
        self.docstring = docstring
        self._default_as_code = None

    def copy_with_replacement(self, name=None, default=_empty, lazy=None,
                              docstring=None):
        if name is None:
            name = self.name
        if default is _empty:
            default = self.default
        if lazy is None:
            lazy = self.lazy
        if docstring is None:
            docstring = self.docstring

        return Parameter(name, default, lazy, docstring)

    def __repr__(self):
        return f"Parameter(name='{self.name}')"

    @property
    def default_as_code(self):
        if self._default_as_code is not None:
            return self._default_as_code
        elif self.default is _empty:
            return None
        elif isinstance(self.default, str):
            # TODO: better escaping here
            return f"'{self.default}'"
        else:
            # this will handle numerics and None
            return str(self.default)

    @default_as_code.setter
    def default_as_code(self, value):
        if not isinstance(value, str):
            raise TypeError("default_as_code must be a string")

        self._default_as_code = value


class Snapshot(StorableObject):
    def __init__(self, **kwargs):
        super().__init__()
        self._reversed_uuid = self.reverse_uuid()  # TODO: change SimStore
        self._reversed_snapshot = None
        for arg, val in kwargs.items():
            setattr(self, arg, val)

    @classmethod
    def _is_compatible(cls, other):
        """check that all required features of ``self`` are in ``other``"""
        required = [p.name for p in cls.__features__.parameters
                    if p.default is p.empty]
        for name in required:
            try:
                _ = getattr(other, name)
            except AttributeError:
                return False

        return True

    def copy_with_replacement(self, **kwargs):
        dct = self.to_dict()

        # copy things that need to be copied
        copy_params = [p for p in self.__features__.parameters
                       if p not in kwargs]
        for param in copy_params:
            dct[param.name] = param.copy(dct[param.name])

        # update with new versions
        dct.update(kwargs)
        return self.from_dict(dct)

    def to_dict(self):
        return {p.name: getattr(self, p.name)
                for p in self.__features__.parameters}

    def _create_reversed(self):
        reversed_dct = {p.name: p.time_reverse(getattr(self, p.name))
                        for p in self.__features__.parameters
                        if 'time_reverse' in p.operations}
        rev = self.copy_with_replacement(**reversed_dct)
        set_uuid(rev, self._reversed_uuid)
        rev._reversed_uuid = get_uuid(self)
        rev._reversed_snapshot = self
        self._reversed_snapshot = rev
        return rev

    @property
    def reverse(self):
        if self._reversed_snapshot is None:
            self._reversed_snapshot = self._create_reversed()
        return self._reversed_snapshot

    @classmethod
    def create_empty(cls):
        new = cls.__new__(cls)
        super(cls, new).__init__()


def minus_reverser(value):
    return -value

def flip_reverser(value):
    return ~value



class function_feature:
    """decorator to label function features"""
    def __init__(self, func):
        self.func = func

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)


class FeatureCollection:
    """

    Note
    ----

        The regular initialization construction should almost never be used
        -- it only creates an empty object. In general, either use the
        ``from_module`` constructor to reate a FeatureCollection from a
        feature module, or combine two feature collections by adding them
        together.
    """
    def __init__(self, parameters=None, functions=None, properties=None,
                 modules=None):
        if parameters is None:
            parameters = []

        if functions is None:
            functions = {}

        if properties is None:
            properties = {}

        if modules is None:
            modules = []

        self.parameters = parameters
        self.functions = functions
        self.properties = properties
        self.modules = modules

    @property
    def classes(self):
        # backwards compatibility TODO: remove
        return self.modules


    @staticmethod
    def _make_parameters_from_old(module):
        parameters = []
        variables = getattr(module, 'variables', [])
        for variable in variables:
            lazy = variable in getattr(module, 'lazy', [])
            if variable in getattr(module, 'default_none', []):
                default = None
            else:
                default = Parameter.empty

            if variable in getattr(module, 'minus', []):
                reverser = minus_reverser
            elif variable in getattr(module, 'flip', []):
                reverser = flip_reverser
            else:
                reverser = None

            if reverser is not None:
                operations = {'time_reverse': reverser}
            else:
                operations = None

            if variable in getattr(module, 'numpy', []):
                copy = np.copy
            else:
                copy = None

            parameters.append(Parameter(
                name=variable,
                default=default,
                lazy=lazy,
                copy=copy,
                operations=operations,
                docstring="TODO",
            ))

        return parameters

    @staticmethod
    def _make_functions_from_old(module):
        functions = getattr(module, 'functions', [])
        return {func: getattr(module, func) for func in functions}


    def __add__(self, other):
        """Combine two FeatureCollections.

        Raise error if there are duplicate entries.
        """
        self._error_if_overlap(other)
        parameters = self.parameters + other.parameters
        functions = {**self.functions, **other.functions}
        properties = {**self.properties, **other.properties}
        modules = self.modules + other.modules
        return self.__class__(parameters, functions, properties, modules)

    @property
    def all_names(self):
        return ([p.name for p in self.parameters] + list(self.properties)
                + list(self.functions))

    def _error_if_overlap(self, other):
        shared = set(self.all_names) & set(other.all_names)
        if shared:
            raise RuntimeError("The following names are used by more than "
                               "one feature for this snapshot: "
                               + str(shared))

    @classmethod
    def from_module(cls, module):

        attrs = ['variables', 'lazy']  # shared module attributes
        # extras will probably not be used in OPS 2.0
        extras = ['required', 'numpy', 'reversal', 'minus', 'flip',
                  'imports', 'debug', 'storables', 'dimensions',
                  'default_none']

        module_dict = vars(module)
        properties = {
            attr: obj for attr, obj in module_dict.items()
            if type(obj) is property
        }

        parameters = [item for item in module_dict.values()
                      if isinstance(item, Parameter)]
        functions = {attr: obj.func for attr, obj in module_dict.items()
                     if isinstance(obj, function_feature)}

        # support building parameters and functions from old-style features
        old_parameters = cls._make_parameters_from_old(module)
        old_functions = cls._make_functions_from_old(module)

        warning = (" is implemented in both old-style and new-style"
                   " snapshot features. Using new style.")

        for parameter in old_parameters:
            if parameter in parameters:
                warnings.warn(f"Parameter named '{parameter.name}'"
                              + warning)
            else:
                parameters.append(parameter)

        for function in old_functions:
            if function in functions:
                warnings.warn(f"Function feature names '{function.name}'"
                              + warning)
            else:
                functions.append(function)

        return cls(parameters, functions, properties, [module])

    def __contains__(self, item):
        return item in self.all_names


def attach_features(features):
    def _decorator(cls):
        # NOTE: this actually updates ``cls`` in place. That's fine when
        # using this with the @-syntax for decorators, but if you use it as
        # a function, that can be a problem:
        # >>> Foo = attach_features(features_1)(Original)
        # >>> Bar = attach_features(features_2)(Original)
        # Bar has features_1 and features_2!

        # gather all the features into one FeatureCollection
        orig = getattr(cls, '__features__', FeatureCollection())
        __features__ = sum([FeatureCollection.from_module(feature)
                            for feature in features], orig)

        # make the docstring
        # TODO

        # methods that we construct (to provide a good signature)
        constructions = {
            '__init__': make_init,
            'copy_with_replacement': make_copy_with_replacement,
        }

        for name, constructor in constructions.items():
            if name not in vars(cls):
                setattr(cls, name, constructor(__features__))
            else:
                pass  # warn that this is user-defined

        # attach properties
        for name, prop in __features__.properties.items():
            setattr(cls, name, prop)

        # attach functions
        for name, func in __features__.functions.items():
            setattr(cls, name, func)

        cls.__features__ = __features__
        return cls

    return _decorator
