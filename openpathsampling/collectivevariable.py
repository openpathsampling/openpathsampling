###############################################################
#| CLASS Order Parameter
###############################################################

import openpathsampling as paths
import chaindict as cd
from openpathsampling.todict import OPSNamed
import marshal
import base64
import types
import opcode

class CollectiveVariable(cd.Wrap, OPSNamed):
    """
    Wrapper for a function that maps a snapshot or an iterable of snapshots
    (like Trajectory) to a number or vector.

    Parameters
    ----------
    name : string
        A descriptive name of the collectivevariable. It is used in the string
        representation.
    dimensions : int
        The number of dimensions of the output order parameter. So far this
        is not used and will be necessary or useful when storage is
        available

    Attributes
    ----------
    name
    dimensions
    single_dict : ChainDict
        The ChainDict that takes care of using only a single element instead of
        an iterable. In the case of a single object. It will be wrapped in a list
        and later only the single element will be returned

    pre_dict : ChainDict
        The ChainDict that will convert all possible input types into parsable
        lists of snapshots, like Trajectory, etc.
    multi_dict
        The ChainDict that will filter duplicate elements for faster calculation
        and compatibility with netcdf format.
    cache_dict
        The ChainDict that will cache calculated values for fast access
    expand_dict
        The ChainDict that will turn unloaded snapshots into real ones so that the
        underlying function will be called with real data
    func_dict
        The ChainDict that will call the actual function in case non of the
        preceeding ChainDicts have returned data

    """

    def __init__(self, name, dimensions=1):
        if (type(name) is not str and type(name) is not unicode) or len(
                name) == 0:
            print type(name), len(name)
            raise ValueError('name must be a non-empty string')

        self.dimensions = dimensions

        self.single_dict = cd.ExpandSingle()
        self.pre_dict = cd.Transform(self._pre_item)
        self.multi_dict = cd.ExpandMulti()
        self.store_dict = cd.MultiStore('collectivevariables', name,
                                        dimensions, self)
        self.cache_dict = cd.ChainDict()

        OPSNamed.__init__(self)

        self.name = name

        if hasattr(self, '_eval'):
            self.expand_dict = cd.UnwrapTuple()
            self.func_dict = cd.Function(None)

            self.func_dict._eval = self._eval

            super(CollectiveVariable, self).__init__(
                post=self.func_dict + self.expand_dict + self.cache_dict +
                     self.store_dict + self.multi_dict + self.single_dict + self.pre_dict
            )

        else:
            super(CollectiveVariable, self).__init__(
                post=self.cache_dict + self.store_dict + self.multi_dict + self.single_dict + self.pre_dict
            )

        self._stored = False

    def flush_cache(self, storage):
        """
        Copy the cache to the internal storage cache for saving

        Parameters
        ----------
        storage : Storage()
            the storage for which the cache should be copied
        """
        if storage in self.store_dict.cod_stores:
            stored = {
                key: value for key, value in self.cache_dict
                if type(key) is tuple or storage in key.idx
            }
            self.store_dict.cod_stores[storage].post.update(stored)
            self.store_dict.cod_stores[storage].update(stored)

    def sync(self, storage):
        """
        Sync this collectivevariable with attached storages

        Parameters
        ----------
        storage : Storage or None
            the store to be used, otherwise all underlying storages are synced
        """
        self.store_dict.update_nod_stores()
        if storage in self.store_dict.cod_stores:
            self.store_dict.cod_stores[storage].sync()

    def cache_all(self, storage):
        """
        Sync this collective variable with attached storages

        Parameters
        ----------
        storage : Storage or None
            the store to be used, otherwise all underlying storages are synced
        """
        self.store_dict.update_nod_stores()
        if storage in self.store_dict.cod_stores:
            self.store_dict.cod_stores[storage].cache_all()

    def _pre_item(self, items):
        item_type = self.store_dict._basetype(items)

        if item_type is paths.Snapshot:
            return items
        elif item_type is paths.Trajectory:
            if len(items) == 0:
                return []
            elif len(items) == 1:
                return [list.__getitem__(items, 0)]
            else:
                return list(list.__iter__(items))
        elif hasattr(items, '__iter__'):
            return list(items)
        else:
            return items

    _compare_keys = ['name', 'dimensions']

    def __eq__(self, other):
        """Override the default Equals behavior"""
        if isinstance(other, self.__class__):
            for key in self._compare_keys:
                if self.__dict__[key] != other.__dict__[key]:
                    return False

            return True

        return NotImplemented

    def __ne__(self, other):
        """Define a non-equality test"""
        if isinstance(other, self.__class__):
            return not self.__eq__(other)
        return NotImplemented

class CV_Volume(CollectiveVariable):
    """ Make `Volume` into `CollectiveVariable`: maps to 0.0 or 1.0

    Attributes
    ----------
    name
    volume
    """

    def __init__(self, name, volume):
        """
        Parameters
        ----------
        name : string
            name of the collective variable
        volume : Volume()
            the Volume instance to be treated as a (storable) CV

        """

        super(CV_Volume, self).__init__(name)
        self.volume = volume

    _compare_keys = ['name', 'volume']

    def _eval(self, items):
        result = [float(self.volume(item)) for item in items]
        return result

    def to_dict(self):
        return {
            'name': self.name,
            'volume': self.volume,
        }

    @classmethod
    def from_dict(cls, dct):
        return CV_Volume(
            name=dct['name'],
            volume=dct['volume']
        )

class CV_Function(CollectiveVariable):
    """Make any function `fcn` into a `CollectiveVariable`.

    Attributes
    ----------
    name
    fcn
    dimensions
    kwargs
    """

    allow_marshal = True
    allowed_modules = ['mdtraj', 'msmbuilder', 'math', 'numpy']

    def __init__(self, name, fcn, dimensions=1, **kwargs):
        """
        Parameters
        ----------
        name : str
        fcn : function
            The function to be used, can almost be an arbitrary function.
        kwargs :
            named arguments which should be given to `fcn` (for example, the
            atoms which define a specific distance/angle). Finally
            `fnc(list_of_snapshots, **kwargs)` is called

        Notes
        -----
        This function is very powerful, but need some explanation if you want
        the function to be stored alongside all other information in your
        storage. The problem is that a python function relies (usually) on
        an underlying global namespace that can be accessed. This is especially
        important when using an iPython notebook. The problem is, that the
        function that stored your used-definited function has no knowledge
        about this underlying namespace and its variables. All it can save is
        names of variables from your namespace to be used. This means you can
        store arbitrary functions, but these will only work, if you call the
        reconstructed ones from the same context (scope). This is a powerful
        feature because a function might do something different in another
        context, but in our case we cannot store these additional information.

        What we can do, is analyse your function and determine which variables
        (if at all these are) and inform you, if you might run into trouble.
        To avoid problem you should try to:
        1. import necessary modules inside of your function
        2. create constants inside your function

        >>> def func(snapshot):
        >>>     import mdtraj as md
        >>>     indices = [4,6,8,10]
        >>>     return md.compute_dihedrals(Trajectory([snapshot]).md(), [indices=indices])

        We will also check if non-standard modules are importet, which are now
        numpy, math, msmbuilder and mdtraj
        """
        super(CV_Function, self).__init__(name, dimensions=dimensions)
        self.callable_fcn = fcn
        self.kwargs = kwargs
        return

    @staticmethod
    def _find_var(code, op):
        opcodes = code.func_code.co_code
        i = 0
        ret = []
        while i < len(opcodes):
            if ord(opcodes[i]) == op:
                ret.append((i, ord(opcodes[i+1]) + ord(opcodes[i+2]) * 256))
            if opcodes[i] < opcode.HAVE_ARGUMENT: #90, as mentioned earlier
                i += 1
            else:
                i += 3

        return [ code.func_code.co_names[i[1]] for i in ret ]

    def to_dict(self):
        fcn = None
        f = self.callable_fcn
        f_module = f.__module__
        is_local = f_module == '__main__'
        is_class = isinstance(f, (type, types.ClassType))
        if not is_class:
            if is_local:
                # this is a local function, let's see if we can save it
                if self.allow_marshal and callable(self.callable_fcn):
                    # use marshal
                    global_vars = CV_Function._find_var(f, opcode.opmap['LOAD_GLOBAL'])
                    import_vars = CV_Function._find_var(f, opcode.opmap['IMPORT_NAME'])

                    if len(global_vars) > 0:
                        print 'Not good. Your function relies on globally set variables ' + \
                              'and these cannot be saved!'
                        print 'requires the following globals to be set:', global_vars
                        print 'Check, if you can replace these by constants or variables ' + \
                              'that are defined within the function itself'

                    not_allowed_modules = [ module for module in import_vars
                                            if module not in CV_Function.allowed_modules]

                    if len(not_allowed_modules) > 0:
                        print 'requires the following modules to be installed:', import_vars
                        print 'note that some of these are not allowed so be careful', not_allowed_modules

                    fcn = {
                        '_marshal': base64.b64encode(
                            marshal.dumps(self.callable_fcn.func_code)),
                        '_global_vars': global_vars,
                        '_module_vars': import_vars
                    }

            elif not is_local:
                # save the external class, e.g. msmbuilder featurizer
                if f_module.split('.')[0] in self.allowed_modules:
                    # only store the function and the module
                    fcn = {
                        '_module': self.callable_fcn.__module__,
                        '_name': self.callable_fcn.__name__
                    }

        return {
            'name': self.name,
            'fcn': fcn,
            'dimensions': self.dimensions,
            'kwargs': self.kwargs
        }

    @classmethod
    def from_dict(cls, dct):
        f = None
        f_dict = dct['fcn']
        if f_dict is not None:
            if '_marshal' in f_dict:
                if cls.allow_marshal:
                    code = marshal.loads(base64.b64decode(f_dict['_marshal']))
                    f = types.FunctionType(code, globals(), code.co_name)

            elif '_module' in f_dict:
                module = f_dict['_module']
                packages = module.split('.')
                if packages[0] in cls.allowed_modules:
                    imp = __import__(module)
                    f = getattr(imp, f_dict['_name'])

        return cls(
            name=dct['name'],
            fcn=f,
            dimensions=dct['dimensions'],
            **dct['kwargs']
        )


    def __eq__(self, other):
        """Override the default Equals behavior"""
        if isinstance(other, self.__class__):
            if self.name != other.name:
                return False
            if self.dimensions != other.dimensions:
                return False
            if self.kwargs != other.kwargs:
                return False

            if hasattr(self.callable_fcn.func_code, 'op_code') and \
                    self.callable_fcn.func_code.op_code != other.callable_fcn.func_code.op_code:
                # Compare Bytecode. Not perfect, but should be good enough
                return False

            return True

        return NotImplemented

    def _eval(self, items):
        trajectory = paths.Trajectory(items)

        return [self.callable_fcn(snap, **self.kwargs) for snap in trajectory]


class CV_Class(CollectiveVariable):
    """Make any callable class `cls` into an `CollectiveVariable`.

    the class instance will be called with single snapshots or a list of
    snapshots. The instance will be created using kwargs.

    Attributes
    ----------
    name
    dimensions
    cls
    kwargs

    Examples
    -------
    >>> # To create an order parameter which calculates the dihedral formed
    >>> # by atoms [7,9,15,17] (psi in Ala dipeptide):
    >>> import mdtraj as md
    >>> traj = None
    >>> psi_atoms = [6,8,14,16]
    >>> psi_orderparam = CV_Function("psi", md.compute_dihedrals,
    >>>                              trajdatafmt="mdtraj",
    >>>                              indices=[psi_atoms])
    >>> print psi_orderparam( traj.md() )
    """

    _allowed_modules = [ 'msmbuilder' ]

    def __init__(self, name, cls, dimensions=1, **kwargs):
        """
        Parameters
        ----------
        name : str
        cls : Class
            a class where instances have a `__call__` attribute
        kwargs : **kwargs
            named arguments which should be given to `cls` (for example, the
            atoms which define a specific distance/angle). Finally an instance
            `instance = cls(**kwargs)` is create when the CV is created and
            using the CV will call `instance(list_of_snapshots)`

        Notes
        -----
        Right now you cannot store user-defined classes. Only classes
        from external packages can be used. For the time being we allow
        classes from msmbuilder.
        """

        super(CV_Class, self).__init__(name, dimensions=dimensions)
        self.callable_cls = cls
        self.kwargs = kwargs
        self._instance = cls(**kwargs)
        return

    def _eval(self, items):
        trajectory = paths.Trajectory(items)

        return [self._instance(snap, **self.kwargs) for snap in trajectory]

    def to_dict(self):
        cls = None
        c = self.callable_cls
        c_module = c.__module__
        is_local = c_module == '__main__'
        is_class = isinstance(c, (type, types.ClassType))

        if is_class:
            if is_local:
                raise TypeError('Locally defined classes are not storable yet')
            elif not is_local:
                if c_module.split('.')[0] in self._allowed_modules:
                    # only store the function and the module
                    cls = {
                        '_module': c.__module__,
                        '_name': c.__name__
                    }

        return {
            'name': self.name,
            'cls': cls,
            'dimensions': self.dimensions,
            'kwargs': self.kwargs
        }

    @classmethod
    def from_dict(cls, dct):
        c = None
        c_dict = dct['fcn']
        if c_dict is not None:
            if '_module' in c_dict:
                module = c_dict['_module']
                packages = module.split('.')
                if packages[0] in cls.allowed_modules:
                    imp = __import__(module)
                    c = getattr(imp, c_dict['_name'])

        return cls(
            name=dct['name'],
            cls=c,
            dimensions=dct['dimensions'],
            **dct['kwargs']
        )

    def __eq__(self, other):
        """Override the default Equals behavior"""
        if isinstance(other, self.__class__):
            if self.name != other.name:
                return False

            if self.fcn != other.fcn:
                return False

            if self.kwargs != other.kwargs:
                return False

            return True

        return NotImplemented



class CV_MD_Function(CV_Function):
    """Make `CollectiveVariable` from `fcn` that takes mdtraj.trajectory as input.

    This is identical to CV_Function except that the function is called with
    an mdraj.Trajetory object instead of the openpathsampling.Trajectory one using
    `fnc(traj.md(), **kwargs)`

    Examples
    -------
    >>> # To create an order parameter which calculates the dihedral formed
    >>> # by atoms [7,9,15,17] (psi in Ala dipeptide):
    >>> import mdtraj as md
    >>> traj = 'paths.Trajectory()'
    >>> psi_atoms = [7,9,15,17]
    >>> psi_orderparam = CV_Function("psi", md.compute_dihedrals,
    >>>                              indices=[[2,4,6,8]])
    >>> print psi_orderparam( traj )
    """

    def __init__(self, name, fcn, dimensions=1, **kwargs):
        """
        Parameters
        ----------
        name : str
        fcn : function
        kwargs :
            named arguments which should be given to `fcn` (for example, the
            atoms which define a specific distance/angle)

        """
        super(CV_MD_Function, self).__init__(
            name, fcn, dimensions=dimensions, **kwargs
        )
        self._topology = None

    def _eval(self, items):
        trajectory = paths.Trajectory(items)

        if self._topology is None:
            self._topology = trajectory.topology.md

        t = trajectory.md(self._topology)
        return self.callable_fcn(t, **self.kwargs)



class CV_Featurizer(CV_Class):
    """
    An CollectiveVariable that uses an MSMBuilder3 featurizer as the logic

    Attribues
    ----------
    name
    cls
    """

    def __init__(self, name, featurizer, **kwargs):
        """

        Parameters
        ----------
        name
        featurizer : msmbuilder.Featurizer
            the featurizer used as a class

        Notes
        -----
        All trajectories or snapshots passed in kwargs will be converted
        to mdtraj objects for convenience
        """

        self.kwargs = kwargs
        # turn Snapshot and Trajectory into md.trajectory
        for key in kwargs:
            if type(kwargs[key]) is paths.Snapshot:
                kwargs[key] = kwargs[key].md()
            elif type(kwargs[key]) is paths.Trajectory:
                kwargs[key] = kwargs[key].md()

        self._feat = featurizer(**kwargs)

        super(CV_Featurizer, self).__init__(
            name,
            cls=featurizer,
            dimensions=self._feat.n_features,
            **self.kwargs
        )


    _compare_keys = ['name']
    _allowed_modules = ['msmbuilder']

    def _eval(self, items):
        trajectory = paths.Trajectory(items)

        # create an MDtraj trajectory out of it
        ptraj = trajectory.md()

        # run the featurizer
        result = self._feat.partial_transform(ptraj)

        return result