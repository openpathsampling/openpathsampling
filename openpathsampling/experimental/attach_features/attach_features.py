from collections import namedtuple
import functools
import copy
import inspect
from numbers import Number

from openpathsampling.netcdfplus import StorableObject

from .write_code import make_init, make_copy_with_replacement

# we use inpect.Parameter.empty in order to more easily convert to
# inspect.Parameter
_empty = inspect.Parameter.empty

class Parameter:
    """Snapshot parameter description.

    Parameters
    ----------
    name: str
        parameter name
    default: Union[Int, Float Str, None]
        default value for this parameter. Must be numeric, string, or None.
        If no default value, use the special value Parameter.empty.
    lazy: bool
        Whether this parameter should use lazy loading
    docstring: str
        The docstring to show for this parmaeter. This should be in numpydoc
        format, and should include the paramter name on the first line.
    """
    empty = _empty
    def __init__(self, name, default=_empty, lazy=False, docstring=None):
        self.name = name
        self.default = default
        self.lazy = lazy
        self.docstring = docstring
        self._default_as_code = None

    def copy_with_replacement(self, name=None, default=None, lazy=None,
                              docstring=None):
        if name is None:
            name = self.name
        if default is None:
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
    def __init__(self):
        super().__init__()
        self._reversed_uuid = None
        self._reversed_snapshot = None
        # TODO: maybe just take kwargs here and have subclass pass kwargs?

    def _is_compatible(self, other):
        # TODO: this should check whether the all the required features of
        # ``other`` are also in ``self``
        pass

    def copy_with_replacement(self, **kwargs):
        dct = self.to_dict()
        dct.update(kwargs)
        return self.from_dict(dct)

    def to_dict(self):
        return {attr: getattr(self, attr)
                for attr in self._snapshot_features.attributes}

    def create_reversed(self):
        rev = self.copy_with_replacement()
        set_uuid(dup, self._reversed_uuid)
        rev._reversed_uuid = get_uuid(self)
        for reverse in self._snapshot_features.reversers:
            reverse(rev)

        rev._reversed = self
        self._reversed = rev
        return rev

    @classmethod
    def create_empty(cls):
        new = cls.__new__(cls)
        super(cls, new).__init__()


def minus_reverser(snap, attr):
    setattr(snap, attr, -getattr(snap, attr))

def flip_reverser(snap, attr):
    setattr(snap, attr, ~getattr(snap, attr))



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
    def __init__(self):
        self.reversers = []  # new addition; func called ``_time_reversal``
        self.parameters = []

        # rename ``classes`` to ``modules``!
        self.classes = []  # used in old storage, new storage

        self.properties = {}
        self.functions = {}  # used in trajectory

    def _make_parameters(self, variables, lazy, default_none):
        # this acts as a bridge between the old approach and the new one
        parameters = []
        for variable in variables:
            lazy = variable in lazy
            if variable in default_none:
                default = None
            else:
                default = Parameter.empty

            parameters.append(Parameter(
                name=variable,
                default=default,
                lazy=lazy,
                # TODO: docstrings?
            ))
        return parameters

    def _make_reversers(self, variables, minus, flip):
        # this acts as a bridge between the old approach and the new one
        reversers = []
        for variable in variables:
            if variable in minus:
                reversers.append(functools.partial(minus_reverser,
                                                   attr=variable))
            if variable in flip:
                reversers.append(functools.partial(flip_reverser,
                                                   attr=variable))
        return reversers

    def __add__(self, other):
        """Combine two FeatureCollections.

        Raise error if there are duplicate entries.
        """
        self._error_if_overlap(other)
        new = copy.deepcopy(self)
        new.parameters += other.parameters
        new.properties.update(other.properties)
        new.functions.update(other.functions)
        return new

    @property
    def modules(self):
        # this will be the name used in the future
        return self.classes

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
        features = cls()

        # new thing
        features.reversers = getattr(module, '_time_reversal', [])

        attrs = ['variables', 'lazy']  # shared module attributes
        # extras will probably not be used in OPS 2.0
        extras = ['required', 'numpy', 'reversal', 'minus', 'flip',
                  'imports', 'debug', 'storables', 'dimensions',
                  'default_none']

        dct = {attr: getattr(module, attr, []) for attr in attrs + extras}


        # identify functions
        features.properties = {
            attr: obj for attr, obj in vars(module).items()
            if type(obj) is property
        }
        func_names = getattr(module, 'functions', [])
        features.functions = {name: getattr(module, name)
                              for name in func_names}

        features.parameters.extend(features._make_parameters(
            dct['variables'], dct['lazy'], dct['default_none']
        ))
        features.reversers.extend(features._make_reversers(
            dct['variables'], dct['flip'], dct['minus']
        ))

        return features

    def __contains__(self, item):
        return item in self.all_names


def _is_function(fnc):
    return callable(fnc) or (hasattr(fnc, "__func__")
                             and callable(fnc.__func__))


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
