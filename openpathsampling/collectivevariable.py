import openpathsampling as paths
import openpathsampling.netcdfplus.chaindict as cd
from openpathsampling.integration_tools import md, error_if_no_mdtraj
from openpathsampling.engines.openmm.tools import trajectory_to_mdtraj
from openpathsampling.netcdfplus import WeakKeyCache, \
    ObjectJSON, create_to_dict, ObjectStore, PseudoAttribute

from openpathsampling.deprecations import (has_deprecations, deprecate,
                                           MSMBUILDER)

import sys
if sys.version_info > (3, ):
    get_code = lambda func: func.__code__
else:
    get_code = lambda func: func.func_code


# ==============================================================================
#  CLASS CollectiveVariable
# ==============================================================================

class CollectiveVariable(PseudoAttribute):
    """
    Wrapper for a function that acts on snapshots or iterables of snapshots

    Parameters
    ----------
    name : string
        A descriptive name of the collectivevariable. It is used in the string
        representation.
    cv_time_reversible : bool
        If ``True`` (default) the CV assumes that reversed snapshots have the
        same value. This is the default case when CVs do not depend on momenta
        reversal. This will speed up computation of CVs by about a factor of
        two. In rare cases you might want to set this to ``False``

    Attributes
    ----------
    name
    cv_time_reversible

    _cache_dict : :class:`openpathsampling.chaindict.ChainDict`
        The ChainDict that will cache calculated values for fast access

    """

    # do not store the settings for the disk cache. These are independent
    # and stored in the cache itself
    _excluded_attr = [
        'diskcache_enabled',
        'diskcache_allow_incomplete',
        'diskcache_chunksize'
    ]

    def __init__(
            self,
            name,
            cv_time_reversible=False
    ):
        super(CollectiveVariable, self).__init__(name, paths.BaseSnapshot)

        self.cv_time_reversible = cv_time_reversible
        self.diskcache_allow_incomplete = not self.cv_time_reversible

        self.diskcache_chunksize = ObjectStore.default_store_chunk_size
        self._cache_dict = cd.ReversibleCacheChainDict(
            WeakKeyCache(),
            reversible=cv_time_reversible
        )

        self._single_dict._post = self._cache_dict

        # self._post = self._single_dict > self._cache_dict

    to_dict = create_to_dict(['name', 'cv_time_reversible'])


class InVolumeCV(CollectiveVariable):
    """Turn a :class:`openpathsampling.volume.Volume` into a collective
    variable

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
        volume : openpathsampling.Volume
            the Volume instance to be treated as a (storable) CV

        """

        super(InVolumeCV, self).__init__(
            name,
            cv_time_reversible=True
        )
        self.volume = volume

        self._eval_dict = cd.Function(
            self._eval,
            requires_lists=False
        )

        self._post = self._post > self._eval_dict

    def _eval(self, items):
        return bool(self.volume(items))

    to_dict = create_to_dict(['name', 'volume'])


class CallableCV(CollectiveVariable):
    """Turn any callable object into a storable :class:`CollectiveVariable`.

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
        cv_callable : callable (function or class with __call__)
            The callable to be used
        cv_time_reversible
        cv_requires_lists : If ``True`` the internal function  always a list of
            elements instead of single values. It also means that if you call
            the CV with a list of snapshots a list of snapshot objects will be
            passed. If ``False`` a list of Snapshots like a trajectory will
            be passed snapshot by snapshot.
        cv_wrap_numpy_array : bool, default: False
            if ``True`` the returned array will be wrapped with a
            ``numpy.array()`` which will convert a list of numpy arrays into a
            single large numpy.array. This is useful for post-processing of
            larger data since numpy arrays are easier to manipulate.
        cv_scalarize_numpy_singletons : bool, default: True
            If ``True`` then arrays of length 1 will be treated as array with
            one dimension less. e.g. [[1], [2], [3]] will be turned into
            [1, 2, 3]. This is often useful, when you use en external function
            to get only a single value.
        **kwargs : kwargs
            a dictionary with named arguments which should be used
            with ``c``. Either for class creation or for calling the function

        Notes
        -----
        This function is abstract and need _eval to be implemented to work.
        Problem is that there are two types of callable functions:
        1. direct functions: these can be called and give the wanted value
           ``c(snapshot, **kwargs)`` would be the typical call
        2. a generating function: a function the creates the callable object
           ``c(**kwargs)(snapshot)`` is the typical call. This is usually used
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
           with the function and this can only be done if they are passed as
           arguments to the function and added as kwargs to the FunctionCV

        >>> import openpathsampling as paths
        >>> def func(snapshot, indices):
        >>>     import mdtraj as md
        >>>     return md.compute_dihedrals(
        >>>         paths.Trajectory([snapshot]).to_mdtraj(), indices=indices)

        >>> cv = FunctionCV('my_cv', func, indices=[[4, 6, 8, 10]])

        The function will also check if non-standard modules are imported,
        which are now ``numpy``, ``math``, ``msmbuilder``, ``pandas`` and
        ``mdtraj``
        """

        super(CallableCV, self).__init__(
            name,
            cv_time_reversible=cv_time_reversible
        )
        self.cv_requires_lists = cv_requires_lists
        self.cv_wrap_numpy_array = cv_wrap_numpy_array
        self.cv_scalarize_numpy_singletons = cv_scalarize_numpy_singletons

        self.cv_callable = cv_callable

        if kwargs is None:
            kwargs = dict()
        self.kwargs = kwargs

        self._eval_dict = cd.Function(
            self._eval,
            self.cv_requires_lists,
            self.cv_scalarize_numpy_singletons
        )

        post = self._post > self._eval_dict

        if cv_wrap_numpy_array:
            # noinspection PyTypeChecker
            post = cd.MergeNumpy() > post

        self._post = post

    def to_dict(self):
        dct = super(CallableCV, self).to_dict()
        callable_argument = self.__class__.args()[2]
        dct[callable_argument] = ObjectJSON.callable_to_dict(self.cv_callable)
        dct['cv_requires_lists'] = self.cv_requires_lists
        dct['cv_wrap_numpy_array'] = self.cv_wrap_numpy_array
        dct['cv_scalarize_numpy_singletons'] = \
            self.cv_scalarize_numpy_singletons
        dct['kwargs'] = self.kwargs
        return dct

    @classmethod
    def from_dict(cls, dct):
        kwargs = dct['kwargs']
        del dct['kwargs']
        dct.update(kwargs)
        obj = cls(**dct)

        return obj

    # def __eq__(self, other):
    #     """Override the default Equals behavior"""
    #     if isinstance(other, self.__class__):
    #         if self.name != other.name:
    #             return False
    #         if self.kwargs != other.kwargs:
    #             return False
    #         if self.cv_callable is None or other.cv_callable is None:
    #             return False
    #
    #         self_code = get_code(self.cv_callable)
    #         other_code = get_code(other.cv_callable)
    #         if hasattr(self_code, 'op_code') \
    #                 and hasattr(other_code, 'op_code') \
    #                 and self_code.op_code != other_code.op_code:
    #             # Compare Bytecode. Not perfect, but should be good enough
    #             return False
    #
    #         return True
    #
    #     return NotImplemented

    __hash__ = CollectiveVariable.__hash__

    def _eval(self, items):
        return items


class FunctionCV(CallableCV):
    """Turn any function into a :class:`CollectiveVariable`.

    Attributes
    ----------
    cv_callable
    """

    def __init__(
            self,
            name,
            f,
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
        cv_time_reversible
        cv_requires_lists
        cv_wrap_numpy_array
        cv_scalarize_numpy_singletons
        **kwargs:
            a dictionary of named arguments which should be given to
            ``cv_callable`` (for example, the atoms which define a specific
            distance/angle). Finally ``cv_callable(snapshots, **kwargs)`` is
            called

        See also
        --------
        :class:`openpathsampling.collectivevariable.CallableCV`

        """

        super(FunctionCV, self).__init__(
            name,
            cv_callable=f,
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


class CoordinateFunctionCV(FunctionCV):
    """Turn any function into a :class:`CollectiveVariable`.

    Attributes
    ----------
    cv_callable
    """

    def __init__(
            self,
            name,
            f,
            cv_requires_lists=False,
            cv_wrap_numpy_array=False,
            cv_scalarize_numpy_singletons=False,
            **kwargs
    ):
        """
        Parameters
        ----------
        name
        f
        cv_requires_lists
        cv_wrap_numpy_array
        cv_scalarize_numpy_singletons
        **kwargs

        See also
        --------
        :class:`openpathsampling.collectivevariable.CallableCV`

        """

        super(FunctionCV, self).__init__(
            name,
            cv_callable=f,
            cv_time_reversible=True,
            cv_requires_lists=cv_requires_lists,
            cv_wrap_numpy_array=cv_wrap_numpy_array,
            cv_scalarize_numpy_singletons=cv_scalarize_numpy_singletons,
            **kwargs
        )

    def to_dict(self):
        dct = super(CoordinateFunctionCV, self).to_dict()
        del dct['cv_time_reversible']
        return dct


class GeneratorCV(CallableCV):
    """Turn a callable class or function generating a callable object into a CV

    The class instance will be called with snapshots. The instance itself
    will be created using the given ``**kwargs``.
    """

    def __init__(
            self,
            name,
            generator,
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
            a class where instances have a ``__call__`` attribute
        cv_time_reversible
        cv_requires_lists
        cv_wrap_numpy_array
        cv_scalarize_numpy_singletons
        **kwargs
            additional arguments which should be given to ``c`` (for example,
            the atoms which define a specific distance/angle). Finally an
            instance ``instance = cls(**kwargs)`` is created when the CV is
            created and using the CV will call ``instance(snapshots)``

        Notes
        -----
        Right now you cannot store user-defined classes. Only classes
        from external packages can be used.
        """

        super(GeneratorCV, self).__init__(
            name,
            cv_callable=generator,
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
        trajectory = paths.Trajectory(items)
        return [self._instance(snap) for snap in trajectory]


class CoordinateGeneratorCV(GeneratorCV):
    """Turn a callable class or function generating a callable object into a CV

    The class instance will be called with snapshots. The instance itself
    will be created using the given ``**kwargs``.
    """

    def __init__(
            self,
            name,
            generator,
            cv_requires_lists=False,
            cv_wrap_numpy_array=False,
            cv_scalarize_numpy_singletons=False,
            **kwargs
    ):
        r"""
        Parameters
        ----------
        name
        generator
        cv_requires_lists
        cv_wrap_numpy_array
        cv_scalarize_numpy_singletons
        **kwargs

        Notes
        -----
        Right now you cannot store user-defined classes. Only classes
        from external packages can be used.
        """

        super(CoordinateGeneratorCV, self).__init__(
            name,
            cv_callable=generator,
            cv_time_reversible=True,
            cv_requires_lists=cv_requires_lists,
            cv_wrap_numpy_array=cv_wrap_numpy_array,
            cv_scalarize_numpy_singletons=cv_scalarize_numpy_singletons,
            **kwargs
        )

    def to_dict(self):
        dct = super(CoordinateGeneratorCV, self).to_dict()
        del dct['cv_time_reversible']
        return dct


class MDTrajFunctionCV(CoordinateFunctionCV):
    """Make ``CollectiveVariable`` from ``f`` that takes
    :class:`mdtraj.Trajectory` as input.

    This is identical to FunctionCV except that the function is called with
    an :class:`mdtraj.Trajectory` object instead of the
    :class:`openpathsampling.Trajectory` one using
    ``f(traj.to_mdtraj(), **kwargs)``

    Examples
    --------
    >>> # To create an order parameter which calculates the dihedral formed
    >>> # by atoms [7,9,15,17] (psi in Ala dipeptide):
    >>> import mdtraj as md
    >>> traj = 'paths.Trajectory()'
    >>> psi_atoms = [7,9,15,17]
    >>> psi_orderparam = FunctionCV("psi", md.compute_dihedrals,
    >>>                              indices=[[2,4,6,8]])
    >>> print psi_orderparam( traj )
    """

    def __init__(self,
                 name,
                 f,
                 topology,
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
        topology : :obj:`openpathsampling.engines.topology.MDTrajTopology`
            the mdtraj topology wrapper from OPS that is used to initialize
            the featurizer in ``pyemma.coordinates.featurizer(topology)``
        cv_requires_lists
        cv_wrap_numpy_array
        cv_scalarize_numpy_singletons
        scalarize_numpy_singletons : bool, default: True
            If ``True`` then arrays of length 1 will be treated as array with
            one dimension less. e.g. ``[[1], [2], [3]]`` will be turned into
            ``[1, 2, 3]``. This is often useful, when you use en external
            function from mdtraj to get only a single value.

        """

        super(MDTrajFunctionCV, self).__init__(
            name,
            f,
            cv_requires_lists=cv_requires_lists,
            cv_wrap_numpy_array=cv_wrap_numpy_array,
            cv_scalarize_numpy_singletons=cv_scalarize_numpy_singletons,
            **kwargs
        )

        self.topology = topology

    def _eval(self, items):
        trajectory = paths.Trajectory(items)

        t = trajectory_to_mdtraj(trajectory, self.topology.mdtraj)
        return self.cv_callable(t, **self.kwargs)

    @property
    def mdtraj_function(self):
        return self.cv_callable

    def to_dict(self):
        return {
            'name': self.name,
            'f': ObjectJSON.callable_to_dict(self.f),
            'topology': self.topology,
            'kwargs': self.kwargs,
            'cv_requires_lists': self.cv_requires_lists,
            'cv_wrap_numpy_array': self.cv_wrap_numpy_array,
            'cv_scalarize_numpy_singletons': self.cv_scalarize_numpy_singletons
        }


@has_deprecations
@deprecate(MSMBUILDER)
class MSMBFeaturizerCV(CoordinateGeneratorCV):
    """A CollectiveVariable that uses an MSMBuilder3 featurizer"""

    def __init__(
            self,
            name,
            featurizer,
            topology,
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
        topology : :obj:`openpathsampling.engines.topology.MDTrajTopology`
            the mdtraj topology wrapper from OPS that is used to initialize
            the featurizer in ``pyemma.coordinates.featurizer(topology)``
        **kwargs :
            a dictionary of named arguments which should be given to ``c``
            (for example, the atoms which define a specific distance/angle).
            Finally an instance ``instance = cls(**kwargs)`` is created when
            the CV is created and using the CV will call
            ``instance(snapshots)``
        cv_wrap_numpy_array
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
            if isinstance(md_kwargs[key], paths.BaseSnapshot):
                md_kwargs[key] = md_kwargs[key].to_mdtraj()
            elif isinstance(md_kwargs[key], paths.Trajectory):
                md_kwargs[key] = md_kwargs[key].to_mdtraj()

        self._instance = featurizer(**md_kwargs)
        self.topology = topology

        super(GeneratorCV, self).__init__(
            name,
            cv_callable=featurizer,
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
        trajectory = paths.Trajectory(items)

        # create an mdtraj trajectory out of it
        ptraj = trajectory_to_mdtraj(trajectory, self.topology.mdtraj)

        # run the featurizer
        return self._instance.partial_transform(ptraj)

    def to_dict(self):
        return {
            'name': self.name,
            'featurizer': ObjectJSON.callable_to_dict(self.featurizer),
            'topology': self.topology,
            'kwargs': self.kwargs,
            'cv_wrap_numpy_array': self.cv_wrap_numpy_array,
            'cv_scalarize_numpy_singletons': self.cv_scalarize_numpy_singletons
        }


class PyEMMAFeaturizerCV(MSMBFeaturizerCV):
    """Make a CV from a function that takes mdtraj.trajectory as input.

    This is identical to :class:`CoordinateGeneratorCV` except that the
    function is called with an :class:`mdraj.Trajetory` object instead of the
    :class:`openpathsampling.Trajectory` one using ``fnc(traj.to_mdtraj(),
    **kwargs)``

    """

    def __init__(
            self,
            name,
            featurizer,
            topology,
            **kwargs
    ):
        """

        Parameters
        ----------
        name
        featurizer : :class:`pyemma.coordinates.featurizer`
            the pyemma featurizer used as a callable class
        topology : :obj:`openpathsampling.engines.topology.MDTrajTopology`
            the mdtraj topology wrapper from OPS that is used to initialize
            the featurizer in ``pyemma.coordinates.featurizer(topology)``
        **kwargs : dict
            a dictionary of named arguments which should be given to the
            ``featurizer`` (for example, the atoms which define a specific
            distance/angle).
            Finally an instance ``instance = cls(**kwargs)`` is created when
            the CV is created and using the CV will call
            ``instance(snapshots)``

        Notes
        -----
        All trajectories or snapshots passed in kwargs will be converted
        to mdtraj objects for convenience
        """

        md_kwargs = dict()
        md_kwargs.update(kwargs)

        # turn Snapshot and Trajectory into md.trajectory
        for key in md_kwargs:
            if isinstance(md_kwargs[key], paths.BaseSnapshot):
                md_kwargs[key] = md_kwargs[key].to_mdtraj()
            elif isinstance(md_kwargs[key], paths.Trajectory):
                md_kwargs[key] = md_kwargs[key].to_mdtraj()

        self.topology = topology

        import pyemma.coordinates
        self._instance = pyemma.coordinates.featurizer(self.topology.mdtraj)

        featurizer(self._instance, **md_kwargs)

        super(GeneratorCV, self).__init__(
            name,
            cv_callable=featurizer,
            cv_requires_lists=True,
            cv_wrap_numpy_array=True,
            cv_scalarize_numpy_singletons=True,
            **kwargs
        )

    def _eval(self, items):
        trajectory = paths.Trajectory(items)

        t = trajectory_to_mdtraj(trajectory, self.topology.mdtraj)
        return self._instance.transform(t)

    def to_dict(self):
        return {
            'name': self.name,
            'featurizer': ObjectJSON.callable_to_dict(self.featurizer),
            'topology': self.topology,
            'kwargs': self.kwargs
        }
