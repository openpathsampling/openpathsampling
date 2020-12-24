from openpathsampling.netcdfplus import DelayedLoader
from .numpydoctools import NumpyDocTools
import openpathsampling as paths

from six import exec_

from collections import namedtuple

import logging

logger = logging.getLogger(__name__)


def _snapshot_function_overridden(cls, method):
    """
    check if in a snapshot class a method was overridden by the user in any (super-)class

    Used to prevent the decorator from overriding user set functions

    Parameters
    ----------
    cls : :class:`BaseSnapshot`
        a class derived from snapshot
    method : str
        the name of the method to be inspected

    Returns
    -------
    bool
        returns `True` if in any super-class the method was overridden

    """

    if hasattr(paths.BaseSnapshot, method):
        fnc = getattr(cls, method)
        # if the method is present in BaseSnapshot it does not count as overridden
        return not (fnc.im_func is getattr(paths.BaseSnapshot, method).im_func)
    elif method in cls.__dict__:
        return True
    elif cls is paths.BaseSnapshot or cls is object:
        return False
    else:
        return _snapshot_function_overridden(cls.__base__, method)


def _register_function(cls, name, code, __features__):

    import numpy as np

    # compile the code and register the new function
    try:
        source_code = '\n'.join(code)
        cc = compile(source_code, '<string>', 'exec')
        #exec cc in locals()
        exec_(cc, locals())

        if name not in cls.__dict__:
            if hasattr(cls, '__features__') and cls.__features__.debug[name] is None:
                raise RuntimeWarning((
                    'Subclassing snapshots with overridden function "%s" is only possible if this '
                    'function is overridden again, otherwise some features might not be copied. '
                    'The general practise of overriding is not recommended.') % name)

            setattr(cls, name, locals()[name])

            __features__['debug'][name] = source_code

        else:
            logger.debug(
                'Function "%s" for class "%s" exists and will not be overridden' %
                (name, cls.__name__)
            )

            __features__['debug'][name] = None

    except RuntimeError as e:
        logger.warn(
            'Problems compiling function "%s" for class "%s". Implementation will not be available!' %
            (name, cls.__name__)
        )
        pass


class CodeFunction(list):
    def __init__(self, name, context):
        list.__init__(self)
        self.name = name
        self.context = context

    def __enter__(self):
        del self[:]
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        _register_function(
            self.context.cls,
            self.name,
            self,
            self.context.__features__
        )

    def __add__(self, other):
        self.extend(other)
        return self

    def format(self, s, list_array_name, include=None, exclude=None):
        if include is None:
            include = []
        if exclude is None:
            exclude = []

        for item in self.context.__features__[list_array_name]:
            add = True
            for inc in include:
                if item not in self.context.__features__[inc]:
                    add = False

            for excl in exclude:
                if item in self.context.__features__[excl]:
                    add = False

            if add:
                self.append(s.format(item))

    def add_uuid(self, name):
        if self.context.use_uuid:
            self += [
                "    {0}.__uuid__ = {0}.get_uuid()".format(name)
            ]

class CodeContext(object):
    def __init__(self, cls, __features__, use_uuid=True):
        self.cls = cls
        self.__features__ = __features__
        self.use_uuid = use_uuid

    def Function(self, name):
        return CodeFunction(name, self)


FeatureTuple = namedtuple(
        'FeatureTuple', 'classes variables properties functions required lazy ' +
                        'numpy reversal minus flip exclude_copy imports debug storables ' +
                        'dimensions default_none'
    )


def attach_features(features, use_lazy_reversed=False):
    """
    Attach features to a snapshot class

    Parameters
    ----------
    features : list of `features`
        a list of features that should be attached to a class. If a class already has
        features attached these will be preserved
    use_lazy_reversed : bool
        if set to `True` the private variables `_reversed` which holds references to
        the reversed snapshot instance will be treated as lazy. This allows the caching
        to get rid of innecessary reversed snapshots to save memory. Otherwise only pairs
        of snapshots can be deleted. Usually you do not need this. But for legacy reasons
        this is still implemented

    """

    # create a parser that can combine numpy docstrings
    parser = NumpyDocTools()

    USE_UUID = True

    def _decorator(cls):
        """
        Class decorator that will attach function for compiled features

        This function will use a list of features and create `__init__` and
        copy functions based on the structure of features for performance.

        It will also take care of creating a joined docstring and the corect
        signature of the `__init__` function

        A attribute `__features__` will be added that contains information
        about the used features their structure. It is a dictionary with the
        following keys

        classes : dict of used features
        lazy : dict of string
            names of features that are treated as lazy loaded object
        properties : dict of string
            names of features that are treated as properties
        lazy : dict of string
            names of features that are treated as lazy loaded object
        reversal : dict of string
            names of features that are treated as being reversible
        minus : dict of string
            names of features that are treated as being reversible and
            should be multiplied by -1.0
        flip : dict of string
            names of features that are treated as being reversible and
            should be negated `~`
        attributes : dict of string
            names of features that are treated as being class attributes
        parameters : dict of string
            names of features that attributes but not properties and hence
            possible parameters for creation
        numpy : dict of string
            names of features that can use numpy for faster copying, etc.

        Parameters
        ----------
        cls : the `class` to the modified

        Returns
        -------
        class
            the modified class

        """

        parser.clear()

        # create and fill `__features__` with values from feature structures
        if hasattr(cls, '__features__'):
            __features__ = {'classes': list(cls.__features__.classes)}
        else:
            __features__ = dict()

        for name in ['variables', 'minus', 'reversal', 'properties',
                     'flip', 'numpy', 'lazy', 'required', 'classes',
                     'exclude_copy', 'imports', 'functions', 'storables',
                     'dimensions', 'default_none']:
            if name not in __features__:
                __features__[name] = []

        if 'debug' not in __features__:
            __features__['debug'] = {}

        for feature in features:
            # check for existing feature and do not register twice
            if feature not in __features__['classes']:
                __features__['classes'].append(feature)

        # add provided additional feature and run though the added ones
        # recursively until nothing is added
        feat_no = 0
        while feat_no < len(__features__['classes']):
            feature = __features__['classes'][feat_no]

            # add provides additional features
            if hasattr(feature, 'imports'):
                if type(feature.imports) is list:
                    for c in feature.imports:
                        if c not in __features__['classes']:
                            __features__['classes'].append(c)
                else:
                    if feature.imports not in __features__['classes']:
                        __features__['classes'].append(feature.imports)

            feat_no += 1

        if use_lazy_reversed:
            cls._reversed = DelayedLoader()

        origin = dict()
        copy_fncs = list()
        copy_feats = list()
        # loop over all the features
        for feature in __features__['classes']:

            # add properties
            for prop in feature.__dict__:
                if hasattr(feature, prop) and type(getattr(feature, prop)) is property:
                    if prop in __features__['properties']:
                        raise RuntimeWarning(
                            'Collision: Property "%s" already exists.' % prop)

                    __features__['properties'] += [prop]
                    setattr(cls, prop, getattr(feature, prop))

            # copy specific attribute types
            for name in ['variables', 'minus', 'lazy', 'flip', 'numpy', 'required', 'imports',
                         'functions', 'storables', 'dimensions', 'default_none']:
                if hasattr(feature, name):
                    content = getattr(feature, name)
                    if type(content) is str:
                        content = [content]

                    def func_or_static(fnc):
                        return callable(fnc) or (hasattr(fnc, "__func__")
                                                 and callable(fnc.__func__))
                    if name == 'functions':
                        for c in content:
                            if hasattr(feature, c):
                                fnc = getattr(feature, c)
                                if func_or_static(fnc):
                                    if hasattr(cls, c):
                                        raise RuntimeWarning(
                                            'Collision: Function "%s" from feature %s already exists.' % (c, feature))
                                    else:
                                        setattr(cls, c, fnc)

                    if name == 'variables':
                        if hasattr(feature, 'copy'):
                            fnc = getattr(feature, 'copy')
                            if callable(fnc):
                                copy_fncs.append(fnc)
                                __features__['exclude_copy'] += content
                                fnc_name = '_copy_' + str(len(copy_feats))
                                setattr(cls, fnc_name, fnc)
                                copy_feats.append(fnc_name)

                        for c in content:
                            if c in __features__['variables']:
                                raise RuntimeError((
                                    'Collision: Attribute "%s" present in two features. ' +
                                    'Please remove one feature "%s" or "%s"') %
                                    (c, str(feature), str(origin[c])))

                    for c in content:
                        origin[c] = feature

                    __features__[name] += content

                else:
                    if name == 'storables':
                        # if storables is missing we assume all variables should be stored

                        __features__['storables'] += __features__['variables']

            # check for cross collisions between variables, properties and function names
            for t1, t2 in [('variables', 'properties'), ('variables', 'functions'), ('properties', 'functions')]:
                union = set(__features__[t1]) & set(__features__[t2])
                if len(union) > 0:
                    raise RuntimeError('Collision: "%s" exist as %s and %s' % (list(union), t1, t2))

        for name in __features__['required']:
            if name not in __features__['variables']:
                raise RuntimeError((
                    'Attribute "%s" is required, but only "%s" are found. Please make sure ' +
                    'that it will be added by ' +
                    'some feature') % (name, str(__features__['variables']))
                )

        # flatten the list of dimensions
        __features__['dimensions'] = list(set(__features__['dimensions']))

        __features__['reversal'] = [
            attr for attr in __features__['variables']
            if attr not in __features__['minus']
            and attr not in __features__['flip']
        ]

        has_lazy = bool(__features__['lazy']) or use_lazy_reversed

        # add descriptors that can handle lazy loaded objects
        for attr in __features__['lazy']:
            setattr(cls, attr, DelayedLoader())

        # update the docstring to be a union of docstrings from the class
        # and the features

        # get docstring from class
        parser.add_docs_from(cls)

        # from top of features
        for feature in __features__['classes']:
            parser.add_docs_from(feature)

            # from properties
            for prop in __features__['properties']:
                if hasattr(feature, prop):
                    if prop not in parser.attributes:
                        parser.add_docs_from(
                            getattr(feature, prop),
                            keep_only=['variables'],
                            translate={'returns': 'variables'}
                        )

        # code for setting default_none (reused in several
        # (for some reason join wasn't working for me?)
        default_none_lines = [
            '    {obj}.' + name + ' = None'
            for name in __features__['default_none']
        ]

        # set new docstring. This is only possible since our class is created
        # using a Metaclass for abstract classes `abc`. Normal classes cannot
        # have their docstring changed.
        cls.__doc__ = parser.get_docstring()

        context = CodeContext(cls, __features__, USE_UUID)

        # compile the function for .copy()

        # def copy(self):
        #     this = cls.__new__(cls)
        #     this._lazy = { ... }
        #     this.feature1 = self.feature1

        with context.Function('copy') as code:
            code += [
                "def copy(self):",
                "    this = cls.__new__(cls)",
            ]

            code.add_uuid('this')
            default_none_code = [line.format(obj='this')
                                 for line in default_none_lines]
            code += default_none_code

            if has_lazy:
                code += [
                    "    this._lazy = {",
                ]
                code.format("       cls.{0} : self._lazy[cls.{0}],",        'lazy', [], ['numpy', 'exclude_copy'])
                code.format("       cls.{0} : self._lazy[cls.{0}].copy(),", 'lazy', ['numpy'], ['exclude_copy'])
                code += [
                    "    }"
                ]

            code += [
                "    this._reversed = None"
            ]

            code.format("    this.{0} = self.{0}",          'variables', [], ['lazy', 'numpy', 'exclude_copy'])
            safe_copy_str  = "    if self.{0} is not None:\n"
            safe_copy_str += "        this.{0} = self.{0}.copy()\n"
            safe_copy_str += "    else:\n"
            safe_copy_str += "        this.{0} = self.{0}"
            # old_copy_str = "    this.{0} = self.{0}.copy()"
            code.format(safe_copy_str,   'variables', ['numpy'], ['lazy', 'exclude_copy'])

            code += map(
                "    self.{0}(this)".format, copy_feats
            )

            code += [
                "    return this"
            ]

        # compile the function for .copy_to(target)

        # def copy_to(self, target):
        #     this = target
        #     this._lazy = { ... }
        #     this.feature1 = self.feature1
        #     return this

        with context.Function('copy_to') as code:
            code += [
                "def copy_to(self, target):",
            ]

            # Copying is effectively creating a new unique object, hence a new UUID
            code.add_uuid('target')

            default_none_code = [line.format(obj='target')
                                 for line in default_none_lines]
            code += default_none_code


            if has_lazy:
                code += [
                    "    target._lazy = {",
                ]
                code.format("       cls.{0} : self._lazy[cls.{0}],",        'lazy', [], ['numpy', 'exclude_copy'])
                code.format("       cls.{0} : self._lazy[cls.{0}].copy(),", 'lazy', ['numpy'], ['exclude_copy'])
                code += [
                    "    }"
                ]

            code += [
                "    target._reversed = None"
            ]

            code.format("    target.{0} = self.{0}",              'variables', [], ['lazy', 'numpy', 'exclude_copy'])
            code.format("    np.copyto(target.{0}, self.{0})",    'variables', ['numpy'], ['lazy', 'exclude_copy'])

            code += map(
                "    self.{0}(this)".format, copy_feats
            )

        # compile the function for .create_reversed()

        # def create_reversed(self):
        #     this = cls.__new__(cls)
        #     this._lazy = { ... }
        #     this.feature1 = self.feature1
        #     this.feature2 = - self.feature2  # minus feature
        #     this.feature3 = ~ self.feature3  # flip features
        #     return this

        with context.Function('create_reversed') as code:
            code += [
                "def create_reversed(self):",
                "    this = cls.__new__(cls)"
            ]

            code += [
                "    this.__uuid__ = self.reverse_uuid()"
            ]

            if has_lazy:
                code += [
                    "    this._lazy = {",
                ]
                code.format("       cls.{0} : self._lazy[cls.{0}],", 'lazy')

                # This should not be necessary since we do not flip lazy loading OPS objects

                # code.format("       cls.{0} : self._lazy[cls.{0}],", 'reversal', ['lazy'])
                # code.format("       cls.{0} : - self._lazy[cls.{0}],", 'minus', ['lazy'])
                # code.format("       cls.{0} : not self._lazy[cls.{0}],", 'flip', ['lazy'])

                code += [
                    "    }"
                ]

            code += [
                "    this._reversed = self"
            ]

            default_none_code = [line.format(obj='this')
                                 for line in default_none_lines]
            code += default_none_code

            code.format("    this.{0} = self.{0}", 'reversal', [], ['lazy'])
            code.format("    this.{0} = - self.{0}", 'minus', [], ['lazy'])
            code.format("    this.{0} = not self.{0}", 'flip', [], ['lazy'])

            code += [
                "    return this"
            ]

        # compile the function for .create_empty()

        # def create_empty(self):
        #     this = cls.__new__(cls)
        #     this._lazy = { ... }
        #     this._reversed = None
        #     return this

        with context.Function('create_empty') as code:
            code += [
                "def create_empty(self):",
                "    this = cls.__new__(cls)"
            ]

            code.add_uuid('this')

            default_none_code = [line.format(obj='this')
                                 for line in default_none_lines]
            code += default_none_code

            if has_lazy:
                code += [
                    "    this._lazy = {}",
                ]

            code += [
                "    this._reversed = None"
            ]
            code += [
                "    return this"
            ]

        # compile the function for __init__

        # def __init__(self, attribute1=None, ... ):
        #     self._lazy = { ... }
        #     self.feature1 = self.feature1
        #     return this

        # we use as signature all feature names in parameters
        parameters = []
        for feature in __features__['variables']:
            if feature in __features__['flip']:
                parameters += ['{0}=False'.format(feature)]
            else:
                parameters += ['{0}=None'.format(feature)]

        if parameters:
            signature = ', ' + ', '.join(parameters)
        else:
            signature = ''

        with context.Function('__init__') as code:
            code += [
                "def __init__(self%s):" % signature,
            ]

            code.add_uuid('self')
            default_none_code = [line.format(obj='self')
                                 for line in default_none_lines]
            code += default_none_code

            # dict for lazy attributes using DelayedLoader descriptor
            if has_lazy:
                code += [
                    "    self._lazy = {",
                ]
                code.format("       cls.{0} : {0},", 'lazy')
                code += [
                    "    }"
                ]

            # set _reversed
            code += [
                "    self._reversed = None"
            ]

            # set non-lazy attributes
            code.format("    self.{0} = {0}", 'variables', [], ['lazy'])

        # compile the function for __init__

        # def init_empty(self)):
        #     self._lazy = {}
        #     self._reversed = None
        #     return this

        with context.Function('init_empty') as code:
            code += [
                "def init_empty(self):",
            ]

            code.add_uuid('self')

            default_none_code = [line.format(obj='self')
                                 for line in default_none_lines]
            code += default_none_code

            # dict for lazy attributes using DelayedLoader descriptor
            if has_lazy:
                code += [
                    "    self._lazy = {}",
                ]

            # set _reversed
            code += [
                "    self._reversed = None"
            ]

        with context.Function('init_copy') as code:
            code += [ "@staticmethod" ]
            code += [
                "def init_copy(self%s):" % signature,
            ]

            code.add_uuid('self')
            default_none_code = [line.format(obj='self')
                                 for line in default_none_lines]
            code += default_none_code

            if has_lazy:
                code += [
                    "    self._lazy = {",
                ]
                code.format("       cls.{0} : {0},",        'lazy', [], ['numpy'])
                code.format("       cls.{0} : {0}.copy(),", 'lazy', ['numpy'], [])
                code += [
                    "    }"
                ]

            code += [
                "    self._reversed = None"
            ]

            code.format("    self.{0} = {0}",          'variables', [], ['lazy', 'numpy'])
            code.format("    np.copyto(self.{0}, {0})",   'variables', ['numpy'], ['lazy'])

        # register (new) __features__ with the class as a namedtuple
        cls.__features__ = FeatureTuple(**__features__)


        return cls

    return _decorator
