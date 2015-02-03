###############################################################
#| CLASS Order Parameter
###############################################################

import mdtraj as md
import openpathsampling as paths
from openpathsampling.todict import restores_as_stub_object

class ObjectDict(dict):
    """
    A cache that is attached to Configuration indices store in the Configuration storage

    Parameters
    ----------
    name : string
        A short and unique name to be used in storage

    Attributes
    ----------
    name

    """

    def __init__(self, dimensions = 1, key_class = None):
        dict.__init__(self)
        self.dimensions = dimensions
        self.key_class = key_class

    def _get(self, obj):
        return dict.__getitem__(self, obj)

    def _set(self, obj, value):
        return dict.__setitem__(self, obj, value)

    def __getitem__(self, items):
        # Allow for numpy style of selecting several indices using a list as index parameter
        if type(items) is list:
            ret = [self._get(key) for key in items]
        else:
            ret = self._get(items)

        return ret

    def __setitem__(self, key, value):
        self._set(key, value)

    def existing(self, objs):
        """
        Find a subset of indices that are present in the cache

        Parameters
        ----------
        indices : list of int
            the initial list of indices to be tested

        Returns
        -------
        existing : list of int
            the subset of indices present in the cache
        """
        return [obj for obj in objs if obj in self]

    def missing(self, objs):
        """
        Find a subset of indices that are NOT present in the cache

        Parameters
        ----------
        indices : list of int
            the initial list of indices to be tested

        Returns
        -------
        existing : list of int
            the subset of indices NOT present in the cache
        """
        return [obj for obj in objs if obj not in self]


class StorableObjectDict(ObjectDict):
    """
    A cache that is attached to Configuration indices store in the Configuration storage

    Parameters
    ----------
    name : string
        A short and unique name to be used in storage

    Attributes
    ----------
    name
    """

    def __init__(self, name, dimensions = 1, key_class = None):
        super(StorableObjectDict, self).__init__(dimensions=dimensions, key_class=key_class)

        self.name = name
        self.storage_caches = dict() # caches the values stored in the data file
        self.var_name = key_class.__name__.lower() + '_' + 'op_' + self.name
        self.object_storages = [] # contains the list of associated ObjectStorages (that itself link to netCDF files)

        for s in self.object_storages:
            if s.content_class is not key_class:
                print 'One of the storages does not store objects of type :', key_class.__name__
                return

        for s in self.object_storages:
            self.storage_caches[s.storage] = dict()

    def storable(self, item):
        """
        Returns True if the given item has indices to be stored in an attached storage otherwise cache it in the dict itself
        """
        for s in self.storage_caches:
            if s in item.idx and item.idx[s] > 0:
                return True
        return False

    def in_store(self, item):
        """
        Returns True, if the item is already stored in an associated cache.
        """
        for s in self.storage_caches:
            if s in item.idx and item.idx[s] in self.storage_caches[s]:
                return True

        return False

    def __contains__(self, item):
        return dict.__contains__(self, item) or self.in_store(item)

    def _get_from_stores(self, item):
        for s in self.storage_caches:
            if s in item.idx and item.idx[s] in self.storage_caches[s]:
                return self.storage_caches[s][item.idx[s]]

        return None

    def _set_to_stores(self, obj, value):
        for s in self.storage_caches:
            if s in obj.idx:
                self.storage_caches[s][obj.idx[s]] = value

    def _set(self, obj, value):
        if self.storable(obj):
            self._set_to_stores(obj, value)
        else:
            return dict.__setitem__(self, obj, value)

    def _get(self, item):
        if self.in_store(item):
            return self._get_from_stores(item)
        else:
            return dict.__getitem__(self, item)

    def __str__(self):
        return "{ 'memory' : " + dict.__str__(self) + ", 'storages' : " + str(self.storage_caches) + " }"


class FunctionalStorableObjectDict(StorableObjectDict):
    """
    A simple dict implementation that will call a function for unknown values

    Parameters
    ----------
    fnc : index (int) -> value (float)
        the function used to generate the cache values for the specified index. In essence a list
    dimensions : int
        the dimension of the stored values. Default is `1`
    content_class : Class
        the class of objects that can be stored
    allow_multiple : boolean
        if True then in the case that the dict is called with several object at a time. The dict
        creates a list of missing ones and passes all of these to the evaluating function at once.
        Otherwise the fall-back is to call each item seperately. If possible always the multiple-
        option should be used.

    Attributes
    ----------
    content_class
    allow_multiple
    dimensions

    """

    def __init__(self, name, fnc, dimensions = 1, key_class = None, allow_multiple = True):
        super(FunctionalStorableObjectDict, self).__init__(
            name=name,
            dimensions=dimensions,
            key_class=key_class
        )
        self._fnc = fnc
        self.allow_multiple = allow_multiple

    def _eval(self, val):
        return self._fnc(val)

    def __call__(self, obj):
        return self._update(obj)

    def _update(self, items):
        if type(items) is list:
            input = items
        else:
            input = [items]

        if self.key_class is not None and len(input) > 0 and isinstance(input[0], self.key_class):
            no_cache = self.missing(input)

            # Add not yet cached data
            if len(no_cache) > 0:
                if self.allow_multiple:
                    result = self._eval(no_cache)

                    for key, res in zip(no_cache, result):
                        self[key] = res
                else:
                    for obj in no_cache:
                        self[obj] = self._eval(obj)

            return self[items]
        else:
            return []

    def get_transformed_view(self, transform):
        def fnc(obj):
            return transform(self(obj))

        return fnc

@restores_as_stub_object
class OrderParameter(FunctionalStorableObjectDict):
    """
    Initializes an OrderParameter object that is essentially a function that maps a frame (Configuration) within a trajectory (Trajectory) to a number.

    Parameters
    ----------
    name : string
        A descriptive name of the orderparameter. It is used in the string representation.
    dimensions : int
        The number of dimensions of the output order parameter. So far this is not used and will be necessary
        or useful when storage is available
    storages : list of ConfigurationStorages()
        contains the list of storages that will be used.

    Attributes
    ----------
    name
    dimensions
    storages

    """

    _instances = dict()

    def __init__(self, name, dimensions = 1):
        if type(name) is str and len(name) == 0:
            raise ValueError('name must be a non-empty string')

#        if name in OrderParameter._instances:
#            raise ValueError(name + ' already exists as an orderparameter. To load an existing one use get_existing(\'' + name + '\')')

        OrderParameter._instances[name] = self
        super(OrderParameter, self).__init__(
            name=name,
            fnc=None,
            dimensions=dimensions,
            key_class=paths.Configuration
        )

    @staticmethod
    def get_existing(name):
        if name in OrderParameter._instances:
            return OrderParameter._instances[name]
        else:
            raise ValueError(name + ' does not exist as an orderparameter')
            return None

    def __eq__(self, other):
        if self is other:
            return True

        if type(other) is type(self):
            if hasattr(self, 'name') and hasattr(other, 'name'):
                if self.name is not None and other.name is not None:
                    return self.name == other.name

        return False

    def __call__(self, items):
        if isinstance(items,  paths.Snapshot):
            return self._update(items.configuration)
        elif isinstance(items, paths.Configuration):
            return self._update(items)
        elif isinstance(items, paths.Trajectory):
            return self._update([snapshot.configuration for snapshot in items])
        elif isinstance(items, list):
            if isinstance(items[0], paths.Configuration):
                return self._update(items)
        else:
            return None

@restores_as_stub_object
class OP_RMSD_To_Lambda(OrderParameter):
    """
    An OrderParameter that transforms the RMSD to a specific center to a lambda value between zero and one.

    Parameters
    ----------
    center : snapshot
        a trajectory snapshot that is used as the point to compute the RMSD to
    lambda_min : float
        rmsd value that corresponds to lambda zero
    max_lambda : float
        rmsd value that corresponds to lambda one
    atom_indices : list of integers (optional)
        a list of integers that is used in the rmsd computation. Usually solvent should be excluded

    Attributes
    ----------
    center : snapshot
        a trajectory snapshot that is used as the point to compute the RMSD to
    lambda_min : float
        rmsd value that corresponds to lambda zero
    max_lambda : float
        rmsd value that corresponds to lambda one
    atom_indices : list of integers (optional)
        a list of integers that is used in the rmsd computation. Usually solvent should be excluded
    metric : msmbuilder.metrics.RMSD
        the RMSD metric object used to compute the RMSD
    _generator : mdtraj.Trajectory prepared by metric.prepare_trajectory
        trajectory object that contains only the center configuration to which the RMSD is computed to
    """

    def __init__(self, name, center, lambda_min, max_lambda, atom_indices=None):
        super(OP_RMSD_To_Lambda, self).__init__(name, dimensions=1)

        self.atom_indices = atom_indices
        self.center = center
        self.min_lambda = lambda_min
        self.max_lambda = max_lambda

        self._generator = paths.Trajectory([center]).subset(self.atom_indices).md()
        return

    ################################################################################
    ##  Actual computation of closest point using RMSD
    ################################################################################

    @staticmethod
    def _scale_fnc(mi, ma):
        def scale(x):
            if x < mi:
                return 0.0
            elif x > ma:
                return 1.0
            else:
                return (x - mi) / (ma - mi)

        return scale

    def _eval(self, items):
        trajectory = paths.Trajectory([paths.Snapshot(configuration=c) for c in items])
        ptraj = trajectory.subset(self.atom_indices).md()

        results = md.rmsd(ptraj, self._generator)

        return map(self._scale_fnc(self.min_lambda, self.max_lambda), results)

@restores_as_stub_object
class OP_Featurizer(OrderParameter):
    """
    An OrderParameter that uses an MSMBuilder3 featurizer as the logic

    Parameters
    ----------
    centers : trajectory
        a trajectory of snapshots that are used as the points to compute the RMSD to
    atom_indices : list of integers (optional)
        a list of integers that is used in the rmsd computation. Usually solvent should be excluded
    metric : msmbuilder.metrics.Metric
        the metric object used to compute the RMSD

    Attributes
    ----------
    centers : trajectory
        a trajectory of snapshots that are used as the points to compute the RMSD to
    atom_indices : list of integers (optional)
        a list of integers that is used in the rmsd computation. Usually solvent should be excluded
    metric : msmbuilder.metrics.RMSD
        the RMSD metric object used to compute the RMSD
    """

    def __init__(self, name, featurizer, atom_indices=None):
        super(OP_Featurizer, self).__init__(name, dimensions=featurizer.n_features)

        self.atom_indices = atom_indices
        self.featurizer = featurizer

        return

    def _eval(self, items):
        trajectory = paths.Trajectory([paths.Snapshot(configuration=c) for c in items])

        # create an MDtraj trajectory out of it
        ptraj = trajectory.subset(self.atom_indices).md()

        # run the featurizer
        result = self.featurizer.partial_transform(ptraj)

        return result

@restores_as_stub_object
class OP_MD_Function(OrderParameter):
    """ Wrapper to decorate any appropriate function as an OrderParameter with a function that need an mdtraj object as input.

    Examples
    -------
    >>> # To create an order parameter which calculates the dihedral formed
    >>> # by atoms [7,9,15,17] (psi in Ala dipeptide):
    >>> import mdtraj as md
    >>> psi_atoms = [7,9,15,17]
    >>> psi_orderparam = OP_Function("psi", md.compute_dihedrals,
    >>>                              indices=[phi_atoms])
    >>> print psi_orderparam( traj.md() )
    """
    def __init__(self, name, fcn, **kwargs):
        """
        Parameters
        ----------
        name : str
        fcn : function
        kwargs :
            named arguments which should be given to `fcn` (for example, the
            atoms which define a specific distance/angle)

        """
        super(OP_MD_Function, self).__init__(name)
        self.fcn = fcn
        self.kwargs = kwargs
        self.topology = None
        return

    def _eval(self, items, *args):
        trajectory = paths.Trajectory(
            [paths.Snapshot(configuration=c) for c in items]
        )

        if self.topology is None:
            # first time ever compute the used topology for this orderparameter to construct the mdtraj objects
            self.topology = trajectory.topology.md

        t = trajectory.md(self.topology)
        return self.fcn(t, *args, **self.kwargs)

@restores_as_stub_object
class OP_Volume(OrderParameter):
    """
    Wrapper that turns a Volume, which can be considered a boolean order parameter into an
    actual OrderParamter()
    """

    def __init__(self, name, volume):
        """
        Parameters
        ----------

        """

        super(OP_Volume, self).__init__(name)
        self.volume = volume

    def _eval(self, items):
        result = [ float(self.volume(item)) for item in items ]
        return result

@restores_as_stub_object
class OP_Function(OrderParameter):
    """ Wrapper to decorate any appropriate function as an OrderParameter.

    Examples
    -------
    >>> # To create an order parameter which calculates the dihedral formed
    >>> # by atoms [7,9,15,17] (psi in Ala dipeptide):
    >>> import mdtraj as md
    >>> psi_atoms = [6,8,14,16]
    >>> psi_orderparam = OP_Function("psi", md.compute_dihedrals,
    >>>                              trajdatafmt="mdtraj",
    >>>                              indices=[psi_atoms])
    >>> print psi_orderparam( traj.md() )
    """
    def __init__(self, name, fcn, **kwargs):
        """
        Parameters
        ----------
        name : str
        fcn : function
        kwargs :
            named arguments which should be given to `fcn` (for example, the
            atoms which define a specific distance/angle)

        Notes
        -----
            We may decide that it is better not to use the trajdatafmt
            trick, and to instead create separate wrapper classes for each
            supported trajformat.
        """
        super(OP_Function, self).__init__(name)
        self._fcn = fcn
        self.kwargs = kwargs
        return


    def _eval(self, items, *args):

        trajectory = paths.Trajectory([paths.Snapshot(configuration=c) for c in items])

        return [ self._fcn(snap, *args, **self.kwargs) for snap in trajectory]
