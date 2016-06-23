###############################################################
# | CLASS Order Parameter
###############################################################

import openpathsampling as paths
import chaindict as cd
from openpathsampling.netcdfplus import StorableNamedObject, WeakLRUCache, ObjectJSON, create_to_dict

import openpathsampling.engines as peng

class CollectiveVariable(cd.Wrap, StorableNamedObject):
    """
    Wrapper for a function that acts on snapshots or iterables of snapshots

    Parameters
    ----------
    name : string
        A descriptive name of the collectivevariable. It is used in the string
        representation.
    cv_time_reversible : bool
        If `True` (default) the CV assumes that reversed snapshots have the same value. This is the
        default case when CVs do not depend on momenta reversal. This will speed up computation of
        CVs by about a factor of two. In rare cases you might want to set this to `False`
    cv_store_cache : bool
        If `True` this CV has a cache on disk attached in form of a table in a netcdf file.
        If set to `False` (default) then there will be no storage created when the cv is stored
        automatically. You can do this later on using function in the cv_store.


    Attributes
    ----------
    name
    cv_time_reversible
    cv_store_cache

    _single_dict : :class:`openpathsampling.chaindict.ChainDict`
        The ChainDict that takes care of using only a single element instead of
        an iterable. In the case of a single object. It will be wrapped in a list
        and later only the single element will be returned
    _cache_dict : :class:`openpathsampling.chaindict.ChainDict`
        The ChainDict that will cache calculated values for fast access

    """

    def __init__(
            self,
            name,
            cv_store_cache=False,
            cv_time_reversible=False
    ):
        if (type(name) is not str and type(name) is not unicode) or len(
                name) == 0:
            raise ValueError('name must be a non-empty string')

        StorableNamedObject.__init__(self)

        self.name = name

        self.cv_store_cache = cv_store_cache
        self.cv_time_reversible = cv_time_reversible

        self._single_dict = cd.ExpandSingle()
        self._cache_dict = cd.ReversibleCacheChainDict(
                WeakLRUCache(1000, weak_type='key'),
                reversible=cv_time_reversible
        )

        super(CollectiveVariable, self).__init__(post=self._single_dict > self._cache_dict)

    def set_cache_store(self, key_store, value_store, backward_store=None):
        """
        Attach store variables to the collective variables.

        If used the collective variable will automatically sync values with
        the store and load from it if necessary. If the CV is created with
        `cv_store_cache = True`. This will be done during CV creation.

        Parameters
        ----------
        key_store : :class:`openpathsampling.netcdfplus.ObjectStore`
            the store that references the key objects used as keys in the
            collective variable (the input for the cv function)
        value_store : :class:`openpathsampling.netcdfplus.ObjectStore`
            the store / variable that references the output (value) objects
        backward_store : :class:`openpathsampling.netcdfplus.ObjectStore` or None
            the optional backward store for reversed objects

        Notes
        -----
        Currently the backward feature is exclusively for
        :class:`openpathsampling.snapshot.BaseSnapshot` which implement a
        reversed object. If `None` backward will equal forward which will
        effectively mean that the store treats forward and backward the
        same.

        """
        self._store_dict = cd.ReversibleStoredDict(
            key_store,
            value_store,
            backward_store,
            self._cache_dict.cache
        )
        self._store_dict._post = self._cache_dict
        self._single_dict._post = self._store_dict

    # This is important since we subclass from list and lists are not hashable
    # but CVs should be
    __hash__ = object.__hash__

    def sync(self):
        """
        Sync this CV with the attached storages

        """
        if hasattr(self, '_store_dict'):
            self._store_dict.sync()

    def cache_all(self):
        """
        Sync this CV with attached storages

        """
        if hasattr(self, '_store_dict'):
            self._store_dict.cache_all()

    def __eq__(self, other):
        """Override the default Equals behavior"""
        if isinstance(other, self.__class__):
            if self.name != other.name:
                    return False

            return True

        return NotImplemented

    def __ne__(self, other):
        """Define a non-equality test"""
        if isinstance(other, self.__class__):
            return not self.__eq__(other)

        return NotImplemented

    to_dict = create_to_dict(['name', 'cv_time_reversible', 'cv_store_cache'])


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
            cv_time_reversible=True
        )
        self.volume = volume

        self._callable_dict = cd.Function(
            self._eval,
            False,
            False
        )

        self._post = self._post > self._callable_dict

    def _eval(self, items):
        result = bool(self.volume(items))
        return result

    to_dict = create_to_dict(['name', 'cv_store_cache', 'volume'])


class CV_Callable(CollectiveVariable):
    """Turn any callable object `cv_callable` into a storable `CollectiveVariable`.

    Attributes
    ----------
    _callable_dict
        The ChainDict that will call the actual function in case non of the
        preceding ChainDicts have returned data
    """

    def __init__(
            self,
            name,
            cv_callable,
            cv_store_cache=True,
            cv_time_reversible=False,
            cv_requires_lists=False,
            cv_wrap_numpy_array=False,
            cv_scalarize_numpy_singletons=False,
            **kwargs
    ):
        r"""
        Parameters
        ----------
        name
        cv_callable : callable (function or class with __call__)
            The callable to be used
        cv_store_cache
        cv_time_reversible
        cv_requires_lists : If `True` the internal function  always a list of elements instead
            of single values. It also means that if you call the CV with a list of snapshots a list
            of snapshot objects will be passed. If `False` a list of Snapshots like a trajectory will
            be passed one by one.
        kwargs : **kwargs
            a dictionary with named arguments which should be used
            with `c`. Either for class creation or for calling the function

        Notes
        -----
        This function is abstract and need _eval to be implemented to work.
        Problem is that there are two types of callable functions:
        1. direct functions: these can be called and give the wanted value
           `c(snapshot, \**kwargs)` would be the typical call
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

        >>> import openpathsampling.engines as peng
        >>> def func(snapshot, indices):
        >>>     import mdtraj as md
        >>>     return md.compute_dihedrals(peng.Trajectory([snapshot]).md(), indices=indices)

        >>> cv = CV_Function('my_cv', func, indices=[[4, 6, 8, 10]])

        The function will also check if non-standard modules are imported, which are now
        numpy, math, msmbuilder, pandas and mdtraj
        """

        super(CV_Callable, self).__init__(
            name,
            cv_store_cache=cv_store_cache,
            cv_time_reversible=cv_time_reversible
        )
        self.cv_requires_lists = cv_requires_lists
        self.cv_wrap_numpy_array = cv_wrap_numpy_array
        self.cv_scalarize_numpy_singletons = cv_scalarize_numpy_singletons

        self.cv_callable = cv_callable

        if kwargs is None:
            kwargs = dict()
        self.kwargs = kwargs

        self._callable_dict = cd.Function(
            self._eval,
            self.cv_requires_lists,
            self.cv_scalarize_numpy_singletons
        )

        post = self._post > self._callable_dict

        if cv_wrap_numpy_array:
            # noinspection PyTypeChecker
            post = cd.MergeNumpy() > post

        self._post = post

    def to_dict(self):
        dct = super(CV_Callable, self).to_dict()
        callable_argument = self.__class__.args()[2]
        dct[callable_argument] = ObjectJSON.callable_to_dict(self.cv_callable)
        dct['cv_requires_lists'] = self.cv_requires_lists
        dct['cv_wrap_numpy_array'] = self.cv_wrap_numpy_array
        dct['cv_scalarize_numpy_singletons'] = self.cv_scalarize_numpy_singletons
        dct['kwargs'] = self.kwargs
        return dct

    @classmethod
    def from_dict(cls, dct):
        kwargs = dct['kwargs']
        del dct['kwargs']
        dct.update(kwargs)
        obj = cls(**dct)

        return obj

    def __eq__(self, other):
        """Override the default Equals behavior"""
        if isinstance(other, self.__class__):
            if self.name != other.name:
                return False
            if self.kwargs != other.kwargs:
                return False
            if self.cv_callable is None or other.cv_callable is None:
                return False

            if hasattr(self.cv_callable.func_code, 'op_code') and hasattr(other.cv_callable.func_code, 'op_code') and \
                            self.cv_callable.func_code.op_code != other.cv_callable.func_code.op_code:
                # Compare Bytecode. Not perfect, but should be good enough
                return False

            return True

        return NotImplemented

    def _eval(self, items):
        return items


class CV_Function(CV_Callable):
    """Turn any function into a `CollectiveVariable`.

    Attributes
    ----------
    cv_callable
    """

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
        r"""
        Parameters
        ----------
        name : str
        f : (callable) function
            The function to be used
        cv_store_cache
        cv_time_reversible
        cv_requires_lists
        cv_wrap_numpy_array
        cv_scalarize_numpy_singletons
        kwargs
            a dictionary of named arguments which should be given to `cv_callable` (for example, the
            atoms which define a specific distance/angle). Finally
            `cv_callable(snapshots, **kwargs)` is called

        See also
        --------
        openpathsampling.CV_Callable

        """

        super(CV_Function, self).__init__(
            name,
            cv_callable=f,
            cv_store_cache=cv_store_cache,
            cv_time_reversible=cv_time_reversible,
            cv_requires_lists=cv_requires_lists,
            cv_wrap_numpy_array=cv_wrap_numpy_array,
            cv_scalarize_numpy_singletons=cv_scalarize_numpy_singletons,
            **kwargs
        )

    @property
    def f(self):
        return self.cv_callable

    def _eval(self, items):
        # here the kwargs are used in the callable when it is evaluated
        return self.cv_callable(items, **self.kwargs)


class CV_Generator(CV_Callable):
    r"""Turn a callable class or other function that generate a callable object into a `CollectiveVariable`.

    The class instance will be called with snapshots. The instance itself
    will be created using the given \**kwargs.
    """

    def __init__(
            self,
            name,
            generator,
            cv_store_cache=True,
            cv_time_reversible=False,
            cv_requires_lists=False,
            cv_wrap_numpy_array=False,
            cv_scalarize_numpy_singletons=False,
            **kwargs
    ):
        r"""
        Parameters
        ----------
        name
        generator : callable class
            a class where instances have a `__call__` attribute
        cv_return_type
        cv_return_shape
        cv_return_simtk_unit
        cv_requires_lists
        cv_store_cache
        kwargs
            additional arguments which should be given to `c` (for example, the
            atoms which define a specific distance/angle). Finally an instance
            `instance = cls(\**kwargs)` is create when the CV is created and
            using the CV will call `instance(snapshots)`

        Notes
        -----
        Right now you cannot store user-defined classes. Only classes
        from external packages can be used.
        """

        super(CV_Generator, self).__init__(
            name,
            cv_callable=generator,
            cv_store_cache=cv_store_cache,
            cv_time_reversible=cv_time_reversible,
            cv_requires_lists=cv_requires_lists,
            cv_wrap_numpy_array=cv_wrap_numpy_array,
            cv_scalarize_numpy_singletons=cv_scalarize_numpy_singletons,
            **kwargs
        )

        # here the kwargs are used when the callable is created (so only once)
        self._instance = generator(**self.kwargs)

    @property
    def instance(self):
        return self._instance

    @property
    def generator(self):
        return self.cv_callable

    def _eval(self, items):
        trajectory = peng.Trajectory(items)
        return [self._instance(snap) for snap in trajectory]


class CV_MDTraj_Function(CV_Function):
    """Make `CollectiveVariable` from `f` that takes mdtraj.trajectory as input.

    This is identical to CV_Function except that the function is called with
    an mdraj.Trajetory object instead of the :class:`openpathsampling.Trajectory` one using
    `f(traj.md(), \**kwargs)`

    Examples
    --------
    >>> # To create an order parameter which calculates the dihedral formed
    >>> # by atoms [7,9,15,17] (psi in Ala dipeptide):
    >>> import mdtraj as md
    >>> traj = 'peng.Trajectory()'
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
        trajectory = peng.Trajectory(items)

        t = trajectory.md()
        return self.cv_callable(t, **self.kwargs)

    @property
    def mdtraj_function(self):
        return self.cv_callable


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
        featurizer : msmbuilder.Featurizer, callable
            the featurizer used as a callable class
        kwargs
            a dictionary of named arguments which should be given to `c` (for example, the
            atoms which define a specific distance/angle). Finally an instance
            `instance = cls(\**kwargs)` is create when the CV is created and
            using the CV will call `instance(snapshots)`
        cv_store_cache
        cv_requires_lists
        cv_scalarize_numpy_singletons

        Notes
        -----
        All trajectories or snapshots passed in kwargs will be converted
        to mdtraj objects for convenience
        """

        md_kwargs = dict()
        md_kwargs.update(kwargs)

        # turn Snapshot and Trajectory into md.trajectory
        for key in md_kwargs:
            if isinstance(md_kwargs[key], peng.BaseSnapshot):
                md_kwargs[key] = md_kwargs[key].md()
            elif isinstance(md_kwargs[key], peng.Trajectory):
                md_kwargs[key] = md_kwargs[key].md()

        self._instance = featurizer(**md_kwargs)

        super(CV_Generator, self).__init__(
            name,
            cv_callable=featurizer,
            cv_store_cache=cv_store_cache,
            cv_time_reversible=True,
            cv_requires_lists=True,
            cv_wrap_numpy_array=cv_wrap_numpy_array,
            cv_scalarize_numpy_singletons=cv_scalarize_numpy_singletons,
            **kwargs
        )

    @property
    def featurizer(self):
        return self.cv_callable

    def _eval(self, items):
        trajectory = peng.Trajectory(items)

        # create an MDtraj trajectory out of it
        ptraj = trajectory.md()

        # run the featurizer
        return self._instance.partial_transform(ptraj)

    def to_dict(self):
        return {
            'name': self.name,
            'featurizer': ObjectJSON.callable_to_dict(self.featurizer),
            'kwargs': self.kwargs,
            'cv_store_cache': self.cv_store_cache,
            'cv_wrap_numpy_array': self.cv_wrap_numpy_array,
            'cv_scalarize_numpy_singletons': self.cv_scalarize_numpy_singletons
        }


class CV_PyEMMA_Featurizer(CV_MSMB_Featurizer):
    """Make `CollectiveVariable` from `fcn` that takes mdtraj.trajectory as input.

    This is identical to CV_Class except that the function is called with
    an mdraj.Trajetory object instead of the openpathsampling.Trajectory one using
    `fnc(traj.md(), **kwargs)`

    """

    def __init__(
            self,
            name,
            featurizer,
            topology,
            cv_store_cache=True,
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
        cv_return_type
        cv_return_shape
        cv_return_simtk_unit
        cv_requires_lists
        cv_store_cache
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
            if isinstance(md_kwargs[key], peng.BaseSnapshot):
                md_kwargs[key] = md_kwargs[key].md()
            elif isinstance(md_kwargs[key], peng.Trajectory):
                md_kwargs[key] = md_kwargs[key].md()

        self.topology = topology

        import pyemma.coordinates
        self._instance = pyemma.coordinates.featurizer(self.topology.md)

        featurizer(self._instance, **md_kwargs)

        super(CV_Generator, self).__init__(
            name,
            cv_callable=featurizer,
            cv_time_reversible=True,
            cv_requires_lists=True,
            cv_store_cache=cv_store_cache,
            cv_wrap_numpy_array=True,
            cv_scalarize_numpy_singletons=True,
            **kwargs
        )

    def _eval(self, items):
        trajectory = peng.Trajectory(items)

        t = trajectory.md(self.topology.md)
        return self._instance.transform(t)

    def to_dict(self):
        return {
            'name': self.name,
            'featurizer': ObjectJSON.callable_to_dict(self.featurizer),
            'topology': self.topology,
            'cv_store_cache': self.cv_store_cache,
            'kwargs': self.kwargs
        }
