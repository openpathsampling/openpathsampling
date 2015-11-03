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
    cv_time_reversible : bool, default: True
        If `True` the CV assumes that reversed snapshots have the same value. This is the
        default case when CVs do not depend on momenta reversal. This will speed up computation of
        CVs by about a factor of two. In rare cases you might want to set this to `False`
    cv_requires_lists : If `True` the internal function  always a list of elements instead
        of single values. It also means that if you call the CV with a list of snapshots a list
        of snapshot objects will be passed. If `False` a list of Snapshots like a trajectory will
        be passed one by one.
    cv_store_cache : bool
        If `True` this CV has a cache on disk attached in form of a table in a netcdf file.
        If set to `False` then there will be no storage created when the cv is stored
        automatically. You can do this later on using function in the cv_store.


    Attributes
    ----------
    name
    cv_time_reversible
    cv_requires_lists
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
            cv_store_cache=False,
            cv_time_reversible=False,
            cv_requires_lists=False,
            cv_wrap_numpy_array=False,
            cv_scalarize_numpy_singletons=False
    ):
        if (type(name) is not str and type(name) is not unicode) or len(
                name) == 0:
            raise ValueError('name must be a non-empty string')

        StorableNamedObject.__init__(self)

        self.name = name

        self.store_cache = cv_store_cache
        self.time_reversible = cv_time_reversible
        self.requires_lists = cv_requires_lists
        self.wrap_numpy_array = cv_wrap_numpy_array
        self.scalarize_numpy_singletons = cv_scalarize_numpy_singletons

        self._single_dict = cd.ExpandSingle()
        self._cache_dict = cd.ReversibleCacheChainDict(WeakLRUCache(1000, weak_type='key'), reversible=cv_time_reversible)

        self._func_dict = cd.Function(
            self._eval,
            self.requires_lists,
            self.scalarize_numpy_singletons
        )

        post = self._func_dict + self._cache_dict + self._single_dict

        if cv_wrap_numpy_array:
            post = post + cd.MergeNumpy()

        super(CollectiveVariable, self).__init__(post=post)

    def set_cache_store(self, key_store, value_store):
        self._store_dict = cd.ReversibleStoredDict(
            key_store,
            value_store,
            self._cache_dict.cache,
            reversible=self.time_reversible
        )
        self._store_dict.post = self._cache_dict
        self._single_dict.post = self._store_dict

    @classmethod
    def from_template(cls, name, f, template, **kwargs):
        f_kwargs = {key: value for key, value in kwargs.iteritems() if not key.startswith('cv_')}
        requires_lists = CollectiveVariable.function_requires_lists(f, template, **f_kwargs)

        cv = cls(name, f, cv_requires_lists=requires_lists, **kwargs)

        # fix cv_returns
        parameters = cls.return_parameters_from_template(f, template, **f_kwargs)
        for key, value in parameters.iteritems():
            setattr(cv, key, value)

        return cv

    # This is important since we subclass from list and lists are not hashable
    # but CVs should be
    __hash__ = object.__hash__

    @staticmethod
    def _identify_var_type(instance):
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

    @staticmethod
    def function_requires_lists(c, template, **kwargs):
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
            `cv_requires_lists` and `cv_return_simtk_unit`

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
            return False

        # Determine if we can use lists or not
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

        return cv_requires_lists

    def return_parameters_from_template(self, template):
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
            `cv_requires_lists` and `cv_return_simtk_unit`

        Notes
        -----
        This is a untility function to create a CV using a template

        See also
        --------
        openpathsampling.CollectiveVariable.from_template
        """

        test_value = self(template)

        cv_return_shape = None
        storable = True
        cv_return_simtk_unit = None

        test_type = test_value

        if type(test_type) is u.Quantity:
            # could be a Quantity([..])
            cv_return_simtk_unit = test_type.unit
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
                cv_return_simtk_unit = test_type.unit
                test_type = test_type._value

        if storable:
            cv_return_type = CollectiveVariable._identify_var_type(test_type)
            return {
                'cv_return_type': cv_return_type,
                'cv_return_shape': cv_return_shape,
                'cv_return_simtk_unit': cv_return_simtk_unit
            }

        return {
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
            'store_cache': self.store_cache,
            'time_reversible': self.time_reversible,
            'requires_lists': self.requires_lists,
            'wrap_numpy_array': self.wrap_numpy_array,
            'scalarize_numpy_singletons': self.scalarize_numpy_singletons
        }

    @classmethod
    def from_dict(cls, dct):
        obj = cls.__new__(cls)
        CollectiveVariable.__init__(
            obj,
            name=dct['name'],
            cv_store_cache=dct['store_cache'],
            cv_time_reversible=dct['time_reversible'],
            cv_requires_lists=dct['requires_lists'],
            cv_wrap_numpy_array=dct['wrap_numpy_array'],
            cv_scalarize_numpy_singletons=dct['scalarize_numpy_singletons']
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
            cv_store_cache=cv_store_cache,
            cv_time_reversible=True,
            cv_requires_lists=True,
            cv_wrap_numpy_array=False,
            cv_scalarize_numpy_singletons=False
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
            cv_store_cache=True,
            cv_time_reversible=False,
            cv_requires_lists=False,
            cv_wrap_numpy_array=False,
            cv_scalarize_numpy_singletons=False,
            **kwargs
    ):
        """
        Parameters
        ----------
        name
        c : callable (function or class with __call__)
            The callable to be used
        cv_store_cache
        cv_time_reversible
        cv_requires_lists
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
            cv_store_cache=cv_store_cache,
            cv_time_reversible=cv_time_reversible,
            cv_requires_lists=cv_requires_lists,
            cv_wrap_numpy_array=cv_wrap_numpy_array,
            cv_scalarize_numpy_singletons=cv_scalarize_numpy_singletons
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
            cv_store_cache=True,
            cv_time_reversible=False,
            cv_requires_lists=False,
            cv_wrap_numpy_array=False,
            cv_scalarize_numpy_singletons=False,
            **kwargs
    ):
        """
        Parameters
        ----------
        name : str
        f : function
            The function to be used
        cv_store_cache
        cv_time_reversible
        cv_requires_lists
        cv_wrap_numpy_array
        cv_scalarize_numpy_singletons
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
            cv_store_cache=cv_store_cache,
            cv_time_reversible=cv_time_reversible,
            cv_requires_lists=cv_requires_lists,
            cv_wrap_numpy_array=cv_wrap_numpy_array,
            cv_scalarize_numpy_singletons=cv_scalarize_numpy_singletons,
            **kwargs
        )

    @property
    def f(self):
        return self.c

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


class CV_Generator(CV_Callable):
    """Turn a callable class or other function that generate a callable object sinto a `CollectiveVariable`.

    The class instance will be called with snapshots. The instance itself
    will be created using the given **kwargs.
    """

    allowed_modules = ['msmbuilder']

    def __init__(
            self,
            name,
            c,
            cv_store_cache=True,
            cv_time_reversible=False,
            cv_requires_lists=False,
            cv_wrap_numpy_array=False,
            cv_scalarize_numpy_singletons=False,
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
        cv_return_simtk_unit
        cv_requires_lists
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

        super(CV_Generator, self).__init__(
            name,
            c=c,
            cv_store_cache=cv_store_cache,
            cv_time_reversible=cv_time_reversible,
            cv_requires_lists=cv_requires_lists,
            cv_wrap_numpy_array=cv_wrap_numpy_array,
            cv_scalarize_numpy_singletons=cv_scalarize_numpy_singletons,
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
        obj = super(CV_Generator, cls).from_dict(dct)
        obj._instance = obj.c(**obj.kwargs)
        return obj


class CV_MDTraj_Function(CV_Function):
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
                 cv_store_cache=True,
                 cv_time_reversible=True,
                 cv_requires_lists=True,
                 cv_wrap_numpy_array=True,
                 cv_scalarize_numpy_singletons=True,
                 **kwargs
                 ):
        """
        Parameters
        ----------
        name : str
        f
        f_kwargs
        cv_store_cache
        cv_time_reversible
        cv_requires_lists
        cv_wrap_numpy_array
        cv_scalarize_numpy_singletons
        scalarize_numpy_singletons : bool, default: True
            If `True` then arrays of length 1 will be treated as array with one dimension less.
            e.g. [ [1], [2], [3] ] will be turned into [1, 2, 3]. This is often useful, when you
            use en external function from mdtraj to get only a single value.

        """

        super(CV_MDTraj_Function, self).__init__(
            name,
            f,
            cv_store_cache=cv_store_cache,
            cv_time_reversible=cv_time_reversible,
            cv_requires_lists=cv_requires_lists,
            cv_wrap_numpy_array=cv_wrap_numpy_array,
            cv_scalarize_numpy_singletons=cv_scalarize_numpy_singletons,
            **kwargs
        )

        self._topology = None

    def _eval(self, items):
        trajectory = paths.Trajectory(items)

        t = trajectory.md()
        return self.f(t, **self.kwargs)


class CV_MSMB_Featurizer(CV_Generator):
    """
    A CollectiveVariable that uses an MSMBuilder3 featurizer

    Attributes
    ----------
    scalarize_numpy_singletons
    """

    def __init__(
            self,
            name,
            featurizer,
            cv_store_cache=True,
            cv_wrap_numpy_array=True,
            cv_scalarize_numpy_singletons=True,
            **kwargs
    ):
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
        cv_store_cache
        cv_requires_lists
        scalarize_numpy_singletons : bool, default: True
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

        super(CV_Generator, self).__init__(
            name,
            c=featurizer,
            cv_store_cache=cv_store_cache,
            cv_time_reversible=True,
            cv_requires_lists=True,
            cv_wrap_numpy_array=cv_wrap_numpy_array,
            cv_scalarize_numpy_singletons=cv_scalarize_numpy_singletons,
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
        return self._instance.partial_transform(ptraj)