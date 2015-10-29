###############################################################
# | CLASS Order Parameter
###############################################################

import marshal
import base64
import types
import opcode
import __builtin__
import importlib

import simtk.unit as u
import numpy as np

import openpathsampling as paths
import chaindict as cd
from openpathsampling.base import StorableNamedObject
from openpathsampling.storage.cache import WeakLRUCache


class CollectiveVariable(cd.Wrap, StorableNamedObject):
    """
    Wrapper for a function that acts on snapshots or iterables of snapshots

    Parameters
    ----------
    name : string
        A descriptive name of the collectivevariable. It is used in the string
        representation.
    cv_return_type : str, default : 'float'
        This specifies the number type of the output of the CV. All types allowed in the netcdfplus.py
        are okay here. Needs to be one of ['bool', 'float', 'index', 'int', 'json', 'lazyobj.*',
        'length', 'long', 'numpy.float32', 'numpy.float64', 'numpy.int16', 'numpy.int32',
        'numpy.int64', 'numpy.int8', 'numpy.uint16', 'numpy.uint32', 'numpy.uint64',
        'numpy.uint8', 'obj.*', 'store', 'str']
    cv_return_shape : None or int or tuple of int, default: None
        A tuple of cv_return_shape of the output of the collective variable.
        `tuple()` corresponds to a scalar, so it `None`. An integer corresponds
        to one dimension of the given length, e.g. `1` corresponds to a one dimensional
        array with length 1. `tuple(1,2,3)` corresponds to a 3-dimensional array of
        size 1 by 2 by 3 elements. The higher dimensional array are usually used with
        numpy arrays.
    cv_requires_lists : If `True` the internal function  always a list of elements instead
        of single values. It also means that if you call the CV with a list of snapshots a list
        of snapshot objects will be passed. If `False` a list of Snapshots like a trajectory will
        be passed one by one.
    cv_simtk_unit : simtk.unit.Unit, default: None
        A simtk.unit.Unit instance specifying the used unit of the output. This means the
        function should return a value with unit. When cached the unit is stripped and when
        loaded recreated.
    cv_store_cache : bool
        If `True` this CV has a cache on disk attached in form of a table in a netcdf file.
        If set to `False` then there will be no storage created when the cv is stored
        automatically. You can do this later on using function in the cv_store.


    Attributes
    ----------
    name
    cv_return_shape
    cv_return_type
    cv_requires_lists
    cv_simtk_unit
    cv_store_cache

    _single_dict : ChainDict
        The ChainDict that takes care of using only a single element instead of
        an iterable. In the case of a single object. It will be wrapped in a list
        and later only the single element will be returned
    _pre_dict : ChainDict
        The ChainDict that will convert all possible input types into
        lists of snapshots, like Trajectory, etc.
    _cache_dict
        The ChainDict that will cache calculated values for fast access
    _func_dict
        The ChainDict that will call the actual function in case non of the
        preceding ChainDicts have returned data

    """

    def __init__(
            self,
            name,
            cv_return_type='float',
            cv_return_shape=None,
            cv_requires_lists=False,
            cv_simtk_unit=None,
            cv_store_cache=False
    ):
        if (type(name) is not str and type(name) is not unicode) or len(
                name) == 0:
            raise ValueError('name must be a non-empty string')

        StorableNamedObject.__init__(self)

        self.name = name

        self.requires_lists = cv_requires_lists
        self.return_shape = cv_return_shape
        self.return_type = cv_return_type
        self.simtk_unit = cv_simtk_unit

        self.store_cache = cv_store_cache

        self._single_dict = cd.ExpandSingle()
        self._cache_dict = cd.CacheChainDict(WeakLRUCache(100000, weak_type='key'))

        self._func_dict = cd.Function(
            self._eval,
            self.requires_lists
        )

        post = self._func_dict + self._cache_dict + self._single_dict

        if 'numpy' in self.return_type:
            post = post + cd.MergeNumpy()

        super(CollectiveVariable, self).__init__(post=post)

    def set_cache_store(self, key_store, value_store):
        self._store_dict = cd.StoredDict(key_store, value_store)
        self._store_dict.post = self._cache_dict
        self._single_dict.post = self._store_dict

    @classmethod
    def from_template(cls, name, f, template, **kwargs):
        f_kwargs = {key: value for key, value in kwargs.iteritems() if not key.startswith('cv_')}
        parameters = cls.parameters_from_template(f, template, **f_kwargs)
        parameters.update(kwargs)
        return cls(name, **parameters)

    # This is important since we subclass from list and lists are not hashable
    # but CVs should be
    __hash__ = object.__hash__

    @staticmethod
    def _interprete_num_type(instance):
        ty = type(instance)

        known_types = [float, int, bool, str]

        if ty in known_types:
            return ty.__name__
        elif hasattr(instance, 'dtype'):
            return 'numpy.' + instance.dtype.type.__name__
        else:
            return 'None'

    def _eval(self, items):
        ### Default CVs don't do anything. Need to use subclass
        return items

    @classmethod
    def parameters_from_template(cls, c, template, **kwargs):
        """
        Compute parameters suitable for a callable using a template snapshot

        Parameters
        ----------
        c : callable (function or class with __call__)
            the callable to be used as a function
        template : openpathsampling.Snapshot
            a test snapshot to be evaluated with the callable

        Returns
        -------
        dict
            A dictionary containing the approriate input parameters for `cv_cv_return_type`, `cv_return_shape`,
            `cv_requires_lists` and `cv_simtk_unit`

        Notes
        -----
        This is a untility function to create a CV using a template

        See also
        --------
        openpathsampling.CollectiveVariable.from_template
        """
        eval_single = True
        value_single = None

        eval_list = True
        value_list = None

        eval_multi = True
        value_multi = None

        cv_requires_lists = None

        try:
            # try use single item
            value_single = c(template, **kwargs)
        except:
            eval_single = False

        try:
            # try use list item
            value_list = c([template], **kwargs)
        except:
            eval_list = False

        try:
            # try use multi list items
            value_multi = c([template, template], **kwargs)
        except:
            eval_multi = False

        if not eval_multi and not eval_list and not eval_single:
            # who knows what happened (after loading), since we
            # cannot use the function we disable the function
            return {
                'c': None
            }

        if eval_single:
            cv_requires_lists = False

        if eval_list and eval_multi:
            # check if results are the same

            #TODO: Check if first and second result are equal. Difficult for numpy
            try:
                if len(value_list) == 1:
                    if len(value_multi) == 2:
                                cv_requires_lists = True

            except TypeError:
                cv_requires_lists = False

        if cv_requires_lists is None:
            # no idea what that function does, but it does not work as
            # expected so we disable it
            return {
                'c': c
            }

        # Determine if storable or not. Means values must be numeric
        # or a list of numeric values

        if cv_requires_lists:
            test_value = value_list[0]
        else:
            test_value = value_single

        cv_return_shape = None
        storable = True
        cv_simtk_unit = None

        test_type = test_value

        if type(test_type) is u.Quantity:
            # could be a Quantity([..])
            cv_simtk_unit = test_type.unit
            test_type = test_type._value

        if type(test_type) is np.ndarray:
            cv_return_shape = test_type.shape
        else:
            if hasattr(test_value, '__len__'):
                cv_return_shape = len(test_value)
                test_type = test_value[0]
                if type(test_type) is u.Quantity:
                    for val in test_value:
                        if type(val._value) is not type(test_value._value):
                            # all values must be of same type
                            storable = False
                else:
                    for val in test_value:
                        if type(val) is not type(test_value):
                            # all values must be of same type
                            storable = False

            if type(test_type) is u.Quantity:
                # could also be [Quantity, ...]
                cv_simtk_unit = test_type.unit
                test_type = test_type._value

        if storable:
            cv_return_type = CollectiveVariable._interprete_num_type(test_type)
            return {
                'c': c,
                'cv_return_type': cv_return_type,
                'cv_return_shape': cv_return_shape,
                'cv_requires_lists': cv_requires_lists,
                'cv_simtk_unit': cv_simtk_unit
            }

        return {
            'c': c
        }

    def sync(self):
        """
        Sync this CV with the attached storages

        Parameters
        ----------
        storage : Storage or None
            the store to be used, otherwise all underlying storages are synced
        """
        if hasattr(self, '_store_dict'):
            self._store_dict.sync()

    def cache_all(self):
        """
        Sync this CV with attached storages

        Parameters
        ----------
        storage : Storage or None
            the store to be used, otherwise all underlying storages are synced
        """
        if hasattr(self, '_store_dict'):
            self._store_dict.cache_all()

    _compare_keys = ['name', 'cv_return_shape']

    def __eq__(self, other):
        """Override the default Equals behavior"""
        if isinstance(other, self.__class__):
            for key in self._compare_keys:
                if getattr(self, key) != getattr(other, key):
                    return False

            return True

        return NotImplemented

    def __ne__(self, other):
        """Define a non-equality test"""
        if isinstance(other, self.__class__):
            return not self.__eq__(other)

        return NotImplemented

    allow_marshal = False
    allowed_modules = []

    @staticmethod
    def _find_var(code, op):
        """
        Helper function to search in python bytecode for a specific function call

        Parameters
        ----------
        code : str
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

    @classmethod
    def callable_to_dict(cls, c):
        """
        Turn a callable function of class into a dictionary

        Used for conversion to JSON

        Parameters
        ----------
        f : callable (function or class with __call__)
            the function to be turned into a dict representaiton
        """
        f_module = c.__module__
        is_local = f_module == '__main__'
        is_loaded = f_module == 'openpathsampling.collectivevariable'
        is_class = isinstance(c, (type, types.ClassType))
        if not is_class:
            if is_local or is_loaded:
                # this is a local function, let's see if we can save it
                if cls.allow_marshal and callable(c):
                    # use marshal
                    global_vars = CV_Function._find_var(c, opcode.opmap['LOAD_GLOBAL'])
                    import_vars = CV_Function._find_var(c, opcode.opmap['IMPORT_NAME'])

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
                                           if module not in CV_Function.allowed_modules]

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
            if f_module.split('.')[0] in cls.allowed_modules:
                # only store the function and the module
                return {
                    '_module': c.__module__,
                    '_name': c.__name__
                }

        raise RuntimeError('Locally defined classes are not storable yet')

    @classmethod
    def callable_from_dict(cls, c_dict):
        """
        Turn a dictionary back in a callable function or class

        Used for conversion from JSON

        Parameters
        ----------
        f_dict : the dictionary that contains the information
        """
        c = None

        if c_dict is not None:
            if '_marshal' in c_dict:
                if cls.allow_marshal:
                    code = marshal.loads(base64.b64decode(c_dict['_marshal']))
                    c = types.FunctionType(code, globals(), code.co_name)

            elif '_module' in c_dict:
                module = c_dict['_module']
                packages = module.split('.')
                if packages[0] in cls.allowed_modules:
                    imp = importlib.import_module(module)
                    c = getattr(imp, c_dict['_name'])

        return c

    def to_dict(self):
        return {
            'name': self.name,
            'return_type': self.return_type,
            'return_shape': self.return_shape,
            'store_cache': self.store_cache,
            'requires_lists': self.requires_lists,
            'simtk_unit': self.simtk_unit
        }

    @classmethod
    def from_dict(cls, dct):
        obj = cls.__new__(cls)
        CollectiveVariable.__init__(
            obj,
            name=dct['name'],
            cv_return_type=dct['return_type'],
            cv_return_shape=dct['return_shape'],
            cv_requires_lists=dct['requires_lists'],
            cv_simtk_unit=dct['simtk_unit'],
            cv_store_cache=dct['store_cache']
        )
        return obj

class CV_Volume(CollectiveVariable):
    """ Turn a `Volume` into a collective variable

    Attributes
    ----------
    name
    volume
    """

    def __init__(self, name, volume, cv_store_cache=True):
        """
        Parameters
        ----------
        name : string
            name of the collective variable
        volume : openpathsampling.Volume
            the Volume instance to be treated as a (storable) CV

        """

        super(CV_Volume, self).__init__(
            name,
            cv_return_type='bool',
            cv_return_shape=None,
            cv_requires_lists=True,
            cv_simtk_unit=None,
            cv_store_cache=cv_store_cache
        )
        self.volume = volume

    _compare_keys = ['name', 'volume']

    def _eval(self, items):
        result = [bool(self.volume(item)) for item in items]
        return result

    def to_dict(self):
        return {
            'name': self.name,
            'volume': self.volume,
            'store_cache': self.store_cache
        }

    @classmethod
    def from_dict(cls, dct):
        return CV_Volume(
            name=dct['name'],
            volume=dct['volume'],
            cv_store_cache=dct['store_cache']
        )


class CV_Callable(CollectiveVariable):
    """Turn any callable object `f` into a storable `CollectiveVariable`.

    Attributes
    ----------
    c
    """

    allow_marshal = True
    allowed_modules = ['mdtraj', 'msmbuilder', 'math', 'numpy', 'pandas']

    def __init__(
            self,
            name,
            c,
            cv_return_type='float',
            cv_return_shape=None,
            cv_requires_lists=False,
            cv_simtk_unit=None,
            cv_store_cache=True,
            **kwargs
    ):
        """
        Parameters
        ----------
        name
        c : callable (function or class with __call__)
            The callable to be used
        cv_return_type
        cv_return_shape
        cv_requires_lists
        cv_simtk_unit
        cv_store_cache
        kwargs : **kwargs
            a dictionary with named arguments which should be used
            with `c`. Either for class creation or for calling the function

        Notes
        -----
        This function is abstract and need _eval to be implemented to work.
        Problem is that there are two types of callable functions:
        1. direct functions: these can be called and give the wanted value
           `c(snapshot, **kwargs)` would be the typical call
        2. a generating function: a function the creates the callable object
           `c(**kwargs)(snapshot)` is the typical call. This is usually used
           for classes. Create the instance and then use it.

        This function is very powerful, but need some explanation if you want
        the function to be stored alongside all other information in your
        storage. The problem is that a python function relies (usually) on
        an underlying global namespace that can be accessed. This is especially
        important when using an iPython notebook. The problem is, that the
        function that stored your used-defined function has no knowledge
        about this underlying namespace and its variables. All it can save is
        names of variables from your namespace to be used. This means you can
        store arbitrary functions, but these will only work, if you call the
        reconstructed ones from the same context (scope). This is a powerful
        feature because a function might do something different in another
        context, but in our case we cannot store these additional information.

        What we can do, is analyse your function and determine which variables
        (if at all these are) and inform you, if you might run into trouble.
        To avoid problems you should try to:
        1. import necessary modules inside of your function
        2. create constants inside your function
        3. if variables from the global scope are used these need to be stored
           with the function and this can only be done if they are passed as arguments
           to the function and added as kwargs to the CV_Function

        >>> def func(snapshot, psi):
        >>>     import mdtraj as md
        >>>     indices = [4,6,8,10]
        >>>     return md.compute_dihedrals(Trajectory([snapshot]).md(), [indices=indices])

        >>> cv = CV_Function('my_cv', func, {'psi': my_global_psi_function})

        The function will also check if non-standard modules are imported, which are now
        numpy, math, msmbuilder, pandas and mdtraj
        """

        super(CV_Callable, self).__init__(
            name,
            cv_return_type=cv_return_type,
            cv_return_shape=cv_return_shape,
            cv_store_cache=cv_store_cache,
            cv_requires_lists=cv_requires_lists,
            cv_simtk_unit=cv_simtk_unit
        )

        self.c = c

        if kwargs is None:
            kwargs = dict()
        self.kwargs = kwargs


    def to_dict(self):
        dct = super(CV_Callable, self).to_dict()
        dct['c'] = self.callable_to_dict(self.c)
        dct['kwargs'] = self.kwargs
        return dct

    @classmethod
    def from_dict(cls, dct):
        obj = super(CV_Callable, cls).from_dict(dct)
        obj.c = cls.callable_from_dict(dct['c'])
        obj.kwargs = dct['kwargs']
        return obj

    def __eq__(self, other):
        """Override the default Equals behavior"""
        if isinstance(other, self.__class__):
            if self.name != other.name:
                return False
            if self.return_shape != other.return_shape:
                return False
            if self.kwargs != other.kwargs:
                return False

            if self.c is None or other.c is None:
                return False

            if hasattr(self.c.func_code, 'op_code') and hasattr(other.c.func_code, 'op_code') and \
                            self.c.func_code.op_code != other.c.func_code.op_code:
                # Compare Bytecode. Not perfect, but should be good enough
                return False

            return True

        return NotImplemented


class CV_Function(CV_Callable):
    """Make any function `f` into a `CollectiveVariable`.

    Attributes
    ----------
    f
    f_kwargs
    """

    allow_marshal = True
    allowed_modules = ['mdtraj', 'msmbuilder', 'math', 'numpy', 'pandas']

    def __init__(
            self,
            name,
            f,
            cv_return_type='float',
            cv_return_shape=None,
            cv_requires_lists=False,
            cv_simtk_unit=None,
            cv_store_cache=True,
            **kwargs
    ):
        """
        Parameters
        ----------
        name : str
        f : function
            The function to be used
        cv_return_type
        cv_return_shape
        cv_requires_lists
        cv_simtk_unit
        cv_store_cache
        kwargs : **kwargs
            a dictionary of named arguments which should be given to `f` (for example, the
            atoms which define a specific distance/angle). Finally
            `f(snapshots, **f_kwargs)` is called

        See also
        --------
        openpathsampling.CV_Callable

        """

        super(CV_Function, self).__init__(
            name,
            c=f,
            cv_return_type=cv_return_type,
            cv_return_shape=cv_return_shape,
            cv_requires_lists=cv_requires_lists,
            cv_simtk_unit=cv_simtk_unit,
            cv_store_cache=cv_store_cache,
            **kwargs
        )

    @property
    def f(self):
        return self.c

    @classmethod
    def parameters_from_template(cls, f, template, **kwargs):
        parameters = super(CV_Function, cls).parameters_from_template(f, template, **kwargs)
        parameters['f'] = parameters['c']
        del parameters['c']
        return parameters

    def _eval(self, items):
        # here the kwargs are used in the callable when it is evaluated
        return self.c(items, **self.kwargs)

    def to_dict(self):
        dct = super(CV_Callable, self).to_dict()
        dct['f'] = self.callable_to_dict(self.f)
        dct['kwargs'] = self.kwargs
        return dct

    @classmethod
    def from_dict(cls, dct):
        obj = super(CV_Callable, cls).from_dict(dct)
        obj.c = cls.callable_from_dict(dct['f'])
        obj.kwargs = dct['kwargs']
        return obj


class CV_Class(CV_Callable):
    """Turn any callable class into a `CollectiveVariable`.

    The class instance will be called with snapshots. The instance itself
    will be created using the given c_kwargs.
    """

    allowed_modules = ['msmbuilder']

    def __init__(
            self,
            name,
            c,
            cv_return_type='float',
            cv_return_shape=None,
            cv_requires_lists=False,
            cv_simtk_unit=None,
            cv_store_cache=True,
            **kwargs
    ):
        """
        Parameters
        ----------
        name
        c : callable class
            a class where instances have a `__call__` attribute
        cv_return_type
        cv_return_shape
        cv_requires_lists
        cv_simtk_unit
        cv_store_cache
        **kwargs : **kwargs
            a dictionary of named arguments which should be given to `c` (for example, the
            atoms which define a specific distance/angle). Finally an instance
            `instance = cls(**kwargs)` is create when the CV is created and
            using the CV will call `instance(snapshots)`

        Notes
        -----
        Right now you cannot store user-defined classes. Only classes
        from external packages can be used.
        """

        super(CV_Class, self).__init__(
            name,
            c=c,
            cv_return_type=cv_return_type,
            cv_return_shape=cv_return_shape,
            cv_store_cache=cv_store_cache,
            cv_requires_lists=cv_requires_lists,
            cv_simtk_unit=cv_simtk_unit,
            **kwargs
        )

        # here the kwargs are used when the callable is created (so only once)
        self._instance = c(**self.kwargs)

    @property
    def instance(self):
        return self._instance

    def _eval(self, items):
        trajectory = paths.Trajectory(items)
        return [self._instance(snap) for snap in trajectory]

    @classmethod
    def from_dict(cls, dct):
        obj = super(CV_Class, cls).from_dict(dct)
        obj._instance = obj.c(**obj.kwargs)
        return obj


class CV_MD_Function(CV_Function):
    """Make `CollectiveVariable` from `f` that takes mdtraj.trajectory as input.

    This is identical to CV_Function except that the function is called with
    an mdraj.Trajetory object instead of the openpathsampling.Trajectory one using
    `f(traj.md(), **kwargs)`

    Examples
    --------
    >>> # To create an order parameter which calculates the dihedral formed
    >>> # by atoms [7,9,15,17] (psi in Ala dipeptide):
    >>> import mdtraj as md
    >>> traj = 'paths.Trajectory()'
    >>> psi_atoms = [7,9,15,17]
    >>> psi_orderparam = CV_Function("psi", md.compute_dihedrals,
    >>>                              indices=[[2,4,6,8]])
    >>> print psi_orderparam( traj )
    """

    def __init__(self,
                 name,
                 f,
                 cv_return_type='numpy.float32',
                 cv_return_shape=None,
                 cv_requires_lists=True,
                 cv_simtk_unit=None,
                 cv_store_cache=True,
                 cv_single_as_scalar=True,
                 **kwargs
                 ):
        """
        Parameters
        ----------
        name : str
        f
        f_kwargs
        cv_return_type
        cv_return_shape
        cv_requires_lists
        cv_simtk_unit
        cv_store_cache
        single_as_scalar : bool, default: True
            If `True` then arrays of length 1 will be treated as array with one dimension less.
            e.g. [ [1], [2], [3] ] will be turned into [1, 2, 3]. This is often useful, when you
            use en external function from mdtraj to get only a single value.

        """

        super(CV_MD_Function, self).__init__(
            name,
            f,
            cv_return_type=cv_return_type,
            cv_return_shape=cv_return_shape,
            cv_requires_lists=cv_requires_lists,
            cv_simtk_unit=cv_simtk_unit,
            cv_store_cache=cv_store_cache,
            **kwargs
        )

        self.single_as_scalar = cv_single_as_scalar
        self._topology = None

    def _eval(self, items):
        trajectory = paths.Trajectory(items)

        if self._topology is None:
            self._topology = trajectory.topology.md

        t = trajectory.md(self._topology)
        arr = self.f(t, **self.kwargs)
        if self.single_as_scalar and arr.shape[-1] == 1:
            return arr.reshape(arr.shape[:-1])
        else:
            return arr

    def to_dict(self):
        dct = super(CV_MD_Function, self).to_dict()
        dct['single_as_scalar'] = self.single_as_scalar
        return dct

    @classmethod
    def from_dict(cls, dct):
        obj = super(CV_MD_Function, cls).from_dict(dct)
        obj.single_as_scalar = dct['single_as_scalar']
        obj._topology = None

        return obj


class CV_MSMB_Featurizer(CV_Class):
    """
    A CollectiveVariable that uses an MSMBuilder3 featurizer

    Attributes
    ----------
    single_as_scalar
    """

    def __init__(self, name, featurizer, cv_store_cache=True, cv_single_as_scalar=True, **kwargs):
        """

        Parameters
        ----------
        name
        c : msmbuilder.Featurizer
            the featurizer used as a callable class
        **kwargs : **kwargs
            a dictionary of named arguments which should be given to `c` (for example, the
            atoms which define a specific distance/angle). Finally an instance
            `instance = cls(**kwargs)` is create when the CV is created and
            using the CV will call `instance(snapshots)`
        cv_return_type
        cv_return_shape
        cv_requires_lists
        cv_simtk_unit
        cv_store_cache
        single_as_scalar : bool, default: True
            If `True` then arrays of length 1 will be treated as array with one dimension less.
            e.g. [ [1], [2], [3] ] will be turned into [1, 2, 3]. This is often useful, when you
            use en external function to get only a single value.

        Notes
        -----
        All trajectories or snapshots passed in kwargs will be converted
        to mdtraj objects for convenience
        """

        md_kwargs = dict()
        md_kwargs.update(kwargs)

        # turn Snapshot and Trajectory into md.trajectory
        for key in md_kwargs:
            if md_kwargs[key].__class__ is paths.Snapshot:
                md_kwargs[key] = md_kwargs[key].md()
            elif md_kwargs[key].__class__ is paths.Trajectory:
                md_kwargs[key] = md_kwargs[key].md()

        self._instance = featurizer(**md_kwargs)
        self.single_as_scalar = cv_single_as_scalar

        return_shape = self._instance.n_features
        if self.single_as_scalar and return_shape == 1:
            return_shape = None

        super(CV_Class, self).__init__(
            name,
            c=featurizer,
            cv_return_shape=return_shape,
            cv_store_cache=cv_store_cache,
            cv_requires_lists=True,
            cv_return_type='numpy.float32',
            cv_simtk_unit=None,
            **kwargs
        )

    @property
    def featurizer(self):
        return self.c

    _compare_keys = ['name']
    allowed_modules = ['msmbuilder']

    def _eval(self, items):
        trajectory = paths.Trajectory(items)

        # create an MDtraj trajectory out of it
        ptraj = trajectory.md()

        # run the featurizer
        arr = self._instance.partial_transform(ptraj)

        if self.single_as_scalar and arr.shape[-1] == 1:
            return arr.reshape(arr.shape[:-1])
        else:
            return arr

    def to_dict(self):
        return {
            'name': self.name,
            'featurizer': self.featurizer,
            'kwargs': self.kwargs,
            'store_cache': self.store_cache,
            'single_as_scalar': self.single_as_scalar
        }

    @classmethod
    def from_dict(cls, dct):
        obj = CV_MSMB_Featurizer(
            name=dct['name'],
            featurizer=dct['featurizer'],
            cv_store_cache=dct['store_cache'],
            cv_single_as_scalar=dct['single_as_scalar'],
            **dct['kwargs']
        )
        return obj