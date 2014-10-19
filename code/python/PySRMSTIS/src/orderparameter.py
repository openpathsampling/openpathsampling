###############################################################
#| CLASS Order Parameter
###############################################################

from msmbuilder.metrics import RMSD
from trajectory import Trajectory
from snapshot import Configuration, Snapshot
import numpy as np

class FunctionalDict(dict):
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

    # TODO: Maybe clean the use of __call__ and __get__ and only use __call__ for the function and
    # then storage in the dict and keep __get__ and __set__ only for regular dict access???

    def __init__(self, fnc, dimensions = 1, content_class = None, allow_multiple = True):
        self._fnc = fnc
        self.content_class = content_class
        self.allow_multiple = allow_multiple
        self.dimensions = dimensions

    def existing(self, items):
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
        return [key for key in items if key in self]

    def missing(self, items):
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
        return [idx for idx in items if idx not in self]

    def _eval(self, val):
        return self._fnc(val)

    def __call__(self, obj):
        return self[obj]

    def _get(self, item):
        return dict.__getitem__(self, item)

    def _update_missing(self, items):
        if type(items) is list:
            input = items
        else:
            input = [items]

        if self.content_class is not None and len(input) > 0 and isinstance(input[0], self.content_class):
            no_cache = self.missing(input)

            # Add not yet cached data
            if len(no_cache) > 0:
                if self.allow_multiple:
                    result = self._eval(no_cache)

                    for key, res in zip(no_cache, result):
                        self[key] = res
                else:
                    for obj in no_cache:
                        self[obj] = self._fnc(obj)
        else:
            return None

    def __getitem__(self, items):
        # Allow for numpy style of selecting several indices using a list as index parameter
        self._update_missing(items)
        if type(items) is list:
            ret = [self._get(key) for key in items]
        else:
            ret = self._get(items)

        return ret


class StorableFunctionDict(FunctionalDict):
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

    def __init__(self, name, fnc, dimensions = 1, content_class = None, storages = None):

        super(StorableFunctionDict, self).__init__(fnc=fnc, dimensions=dimensions, content_class=content_class)

        self.name = name
        self.storage_caches = dict() # caches the values stored in the data file
        self.var_name = content_class.__name__.lower() + '_' + 'op_' + self.name
        self.object_storages = [] # contains the list of associated ObjectStorages (that itself link to netCDF files)

        if storages is None:
            self.object_storages = []
        elif type(storages) is list:
            self.object_storages = storages
        else:
            self.object_storages = [storages]

        for s in self.object_storages:
            if s.content_class is not content_class:
                print 'One of the storages does not store objects of type :', content_class.__name__
                return

        for s in self.object_storages:
            self.storage_caches[s.storage] = dict()

        self._init_storage()
        self._update_store()


    def in_store(self, item):
        """
        Returns True, if the item is already stored in an associated cache.
        """
        for s in self.storage_caches:
            if s in item.idx and item.idx[s] in self.storage_caches[s]:
                return True

        return False

    def __contains__(self, item):
        return super(StorableFunctionDict, self).__contains__(item) or self.in_store(item)

    def _get_from_stores(self, item):
        for s in self.storage_caches:
            if s in item.idx and item.idx[s] in self.storage_caches[s]:
                return self.storage_caches[s][item.idx[s]]

        return None

    def _get(self, item):
        if self.in_store(item):
            return self._get_from_stores(item)
        else:
            return dict.__getitem__(self, item)

    def _init_storage(self, object_storage = None):
        """
        initializes the associated storage to index a specific order parameter in it
        """
        if object_storage is None:
            # just run this function with all registered storages
            if len(self.object_storages) > 0:
                map(self._init_storage, self.object_storages)
        else:
            storage = object_storage.storage
            if self.var_name not in storage.variables:
                # define dimensions for non-scalar orderparameters
                size_dimension_name = 'scalar'
                if self.dimensions > 1:
                    size_dimension_name = self.var_name + "_dimension"
                    storage.createDimension(size_dimension_name, self.dimensions)

                idx_dimension = self.object_storages[0].idx_dimension

                # define variables for OrderParameters
                ncvar_op = storage.createVariable(self.var_name, 'f', (idx_dimension, size_dimension_name))
                ncvar_op_idx = storage.createVariable(self.var_name + '_idx', 'u4', idx_dimension)
                ncvar_op_length = storage.createVariable(self.var_name + '_length', 'u4', 'scalar')

                ncvar_op_length[0] = 0

                # Define units for configuration variables.
                setattr(ncvar_op, 'units', 'None')

                # Define long (human-readable) names for variables.
                setattr(ncvar_op, "long_name", "op_" + self.name + "[" + idx_dimension + "][" + size_dimension_name + "] is the orderparameter '" + self.name + "' of configuration 'configuration'.")
            else:
                self.load()

    def _update_store(self, storage = None):
        """
        This will transfer everything from the memory cache into the storage copy in memory which is used to interact with
        the file storage.

        Parameters
        ----------
        storage : Storage() on None
            The storage (not ObjectStorage) to store in. If None then all associated storages will be updated up.

        """
        if storage is None:
            if len(self.storage_caches) > 0:
                map(self._update_store, self.storage_caches.keys())
        else:
            if storage not in self.storage_caches:
                # TODO: Throw exception
                self.storage_caches[storage] = dict()

            store = self.storage_caches[storage]
            for item, value in self.iteritems():
                if storage in item.idx:
                    store[item.idx[storage]] = value

    def tidy_cache(self, storage = None):
        """
        This will transfer everything from the memory cache into the storage copy in memory which is used to interact with
        the file storage.

        Parameters
        ----------
        storage : Storage() on None
            The storage (not ObjectStorage) to store in. If None then all associated storages will be cleaned up.

        """
        if storage is None:
            if len(self.storage_caches) > 0:
                map(self.tidy_cache, self.storage_caches.keys())
        else:
            # Make sure configuration_indices are stored and have an index and then _add_class the configuration index to the trajectory

            if storage not in self.storage_caches:
                # TODO: Throw exception
                self.storage_caches[storage] = dict()

            for item, value in self.iteritems():
                if storage in item.idx:
                    del self[item]

    def save(self, storage = None):
        """
        Save the current state of the cache to the storage.

        Parameters
        ----------
        storage : Storage() on None
            The storage (not ObjectStorage) to store in. If None then all associated storages will be saved in.

        """

        if storage is None:
            if len(self.storage_caches) > 0:
                map(self.save, self.storage_caches.keys())
                map(self.tidy_cache, self.storage_caches.keys())
        else:
            # Make sure configuration_indices are stored and have an index and then _add_class the configuration index to the trajectory

            self._update_store(storage)
            store = self.storage_caches[storage]

            storage.variables[self.var_name][:, :] = np.array(store.values())
            storage.variables[self.var_name + '_idx'][:] = np.array(store.keys())
            storage.variables[self.var_name + '_length'][0] = len(store)

    def load(self, storage = None):
        """
        Restores the cache from the storage using the name of the orderparameter.

        Parameters
        ----------
        storage : Storage() on None
            The storage (not ObjectStorage) to store in. If None then all associated storages will be loaded from.

        Notes
        -----
        Make sure that you use unique names otherwise you might load the wrong parameters!
        """

        if storage is None:
            if len(self.storage_caches) > 0:
                map(self.load, self.storage_caches.keys())
        else:
            length = int(storage.variables[self.var_name + '_length'][0])
            data_idx = storage.variables[self.var_name + '_idx'][:length].astype(np.int).copy()
            data = storage.variables[self.var_name][:, 0].astype(np.float).tolist()
            self.storage_caches[storage] = dict(zip(data_idx, data))


class OrderParameter(StorableFunctionDict):
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

    def __init__(self, name, dimensions = 1, storages = None):
        super(OrderParameter, self).__init__(name, None, dimensions, Configuration, storages=storages)

    def __call__(self, items):
        if isinstance(items, Snapshot):
            return self[items.configuration]
        elif isinstance(items, Configuration):
            return self[items]
        elif isinstance(items, Trajectory):
            return self[[snapshot.configuration for snapshot in items]]
        elif isinstance(items, list):
            if isinstance(items[0], Configuration):
                return self[items]

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

    def __init__(self, name, center, lambda_min, max_lambda, atom_indices=None, storages = None):
        super(OP_RMSD_To_Lambda, self).__init__(name, dimensions=1, storages=storages)

        self.atom_indices = atom_indices
        self.center = center
        self.min_lambda = lambda_min
        self.max_lambda = max_lambda

        # Generate RMSD metric using only the needed indices. To index memory we crop the read snapshots from the database and do not use a cropping RMSD on the full snapshots
        self.metric = RMSD(None)
        self._generator = self.metric.prepare_trajectory(Trajectory([center]).subset(self.atom_indices).md())
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
        trajectory = Trajectory([Snapshot(configuration=c) for c in items])

        ptraj = self.metric.prepare_trajectory(trajectory.subset(self.atom_indices).md())
        results = self.metric.one_to_all(self._generator, ptraj, 0)

        return map(self._scale_fnc(self.min_lambda, self.max_lambda), results)


class OP_Multi_RMSD(OrderParameter):
    """
    An OrderParameter that transforms the RMSD to a set of centers.

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

    def __init__(self, name, centers, atom_indices=None, metric=None, storages=None):
        super(OP_Multi_RMSD, self).__init__(name, dimensions=len(centers), storages=storages)

        self.atom_indices = atom_indices
        self.center = centers

        # Generate RMSD metric using only the needed indices. To index memory we crop the read snapshots from the database and do not use a cropping RMSD on the full snapshots
        if metric is None:
            self.metric = RMSD(None)
        else:
            self.metric = metric
        self._generator = self.metric.prepare_trajectory(centers.subset(self.atom_indices).md())
        return

    def _eval(self, items):
        trajectory = Trajectory([Snapshot(configuration=c) for c in items])

        ptraj = self.metric.prepare_trajectory(trajectory.subset(self.atom_indices).md())
        return [self.metric.one_to_all(ptraj, self._generator, idx) for idx in range(0, len(ptraj))]

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
        super(OP_Function, self).__init__(name)
        self.fcn = fcn
        self.kwargs = kwargs
        self.topology = None
        return

    def _eval(self, trajectory, *args):
        if self.topology is None:
            # first time ever compute the used topology for this orderparameter to construct the mdtraj objects
            self.topology = trajectory.md_topology()

        t = trajectory.md(self.topology)
        return self.fcn(t, *args, **self.kwargs)


class OP_Function(OrderParameter):
    """ Wrapper to decorate any appropriate function as an OrderParameter.

    Examples
    -------
    >>> # To create an order parameter which calculates the dihedral formed
    >>> # by atoms [7,9,15,17] (psi in Ala dipeptide):
    >>> import mdtraj as md
    >>> psi_atoms = [7,9,15,17]
    >>> psi_orderparam = OP_Function("psi", md.compute_dihedrals,
    >>>                              trajdatafmt="mdtraj",
    >>>                              indices=[phi_atoms])
    >>> print psi_orderparam( traj.md() )
    """
    def __init__(self, name, fcn, trajdatafmt=None, **kwargs):
        """
        Parameters
        ----------
        name : str
        fcn : function
        trajdatafmt : str
            which format the trajectory data needs to be in for the `fcn`.
            Currently supports "mdtraj", otherwise defaults to our own
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
        self.fcn = fcn
        self.trajdatafmt = trajdatafmt
        self.kwargs = kwargs
        return


    def _eval(self, trajectory, *args):
        if self.trajdatafmt=='mdtraj':
            t = trajectory.md()
        else:
            t = trajectory
        return self.fcn(t, *args, **self.kwargs)


if __name__ == '__main__':
    def ident(indices):
        if type(indices) is list:
            return [float(i) for i in indices]
        else:
            return float(indices)

    s = StorableFunctionDict(name='TestList', fnc=ident)

    print s[10]
    print s[5]
    print s[5]
    print s.cache
    print s[1:10]
    print s.cache