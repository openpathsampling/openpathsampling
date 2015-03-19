###############################################################
#| CLASS Order Parameter
###############################################################

import mdtraj as md
import openpathsampling as paths
from openpathsampling.todict import restores_as_stub_object

class NestableObjectDict(dict):
    """
    Cache attached to Configuration indices stored in Configuration storage

    Parameters
    ----------
    name : string
        A short and unique name to be used in storage

    Attributes
    ----------
    name

    fnc : index (int) -> value (float)
        the function used to generate the cache values for the specified index. In essence a list
    dimensions : int
        the dimension of the stored values. Default is `1`
    content_class : Class
        the class of objects that can be stored
    fnc_uses_lists : boolean
        if True then in the case that the dict is called with several object at a time. The dict
        creates a list of missing ones and passes all of these to the evaluating function at once.
        Otherwise the fall-back is to call each item seperately. If possible always the multiple-
        option should be used.

    Attributes
    ----------
    content_class
    fnc_uses_lists
    dimensions
    """

    use_unique = True

    def __init__(self):
        dict.__init__(self)
        self.post = None

    def __iter__(self):
        return None

    def __getitem__(self, items):
        is_list = type(items) is list

        if is_list:
            results = self._get_list(items)
        else:
            results = self._get(items)

        if self.post is not None:
            if is_list:
                nones = [obj[0] for obj in zip(items, results) if obj[1] is None]
                if len(nones) == 0:
                    return results
                else:
                    rep = self.post[[p for p in nones]]
                    self._add_new(nones, rep)

                    it = iter(rep)
                    return [it.next() if p[1] is None else p[1] for p in zip(items, results)]
            else:
                if results is None:
                    rep = self.post[items]

                    self._add_new(items, rep)
                    return rep
                else:
                    return results

        return results

    def _add_new(self, items, values):
        self[items] = values

    def __setitem__(self, key, value):
        if type(key) is list:
            self._set_list(key, value)
        else:
            self._set(key, value)

#    def __contains__(self, item):
#        return dict.__contains__(self, item) or self.in_store(item)


    def _contains(self, item):
        return dict.__contains__(self, item)

    def _contains_list(self, items):
        return [dict.__contains__(self, item) for item in items]

    def _set(self, item, value):
        dict.__setitem__(self, item, value)

    def _set_list(self, items, values):
        [dict.__setitem__(self, item, value) for item, value in zip(items, values)]

    def _get(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            return None

    def _get_list(self, items):
        return [ self._get(item) for item in items ]

    def push(self):
        pass


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

    def __call__(self, items):
        return self[items]

    def __add__(self, other):
        other.post = self
        return other

    def _split_list_dict(self, dct, items):
        nones = [item[self] if item in self else None for item in items]
        missing = [item for item in items if item in self]

        return nones, missing

    def _split_list(self, keys, values):
        missing = [ obj[0] for obj in zip(keys,values) if obj[1] is None ]
        nones = values

        return nones, missing

    def _apply_some_list(self, func, items):
        some = [item for item in items if item is not None]
        replace = func(some)
        it = iter(replace)
        return [ it.next() if obj is not None else None for obj in some ]

    def _replace_none(self, nones, replace):
        it = iter(replace)
        return [ obj if obj is not None else it.next() for obj in nones ]


class NODExpandMulti(NestableObjectDict):
    """
    Will only request the unique keys to post
    """

    def __getitem__(self, items):
        is_list = type(items) is list

        if not is_list:
            return self.post._get(items)

        if len(items) == 0:
            return []

        uniques = list(set(items))
        rep_unique = self.post[[p for p in uniques]]
        multi_cache = dict(zip(uniques, rep_unique))

        return [multi_cache[item] for item in items]

    def __setitem__(self, key, value):
        self.post[key] = value

    def _add_new(self, items, values):
        pass

class NODTransform(NestableObjectDict):
    def __init__(self, transform, post):
        super(NODTransform, self).__init__()
        self.transform = transform

    def __getitem__(self, item):
        return self.post[self.transform(item)]

    def __setitem__(self, key, value):
        self.post[self.transform(key)] = value

    def _add_new(self, items, values):
        pass


class NODFunction(NestableObjectDict):
    def __init__(self, fnc, fnc_uses_lists=True):
        super(NODFunction, self).__init__()
        self._fnc = fnc
        self.fnc_uses_lists = fnc_uses_lists

    def _eval(self, val):
        return self._fnc(val)

    def _contains(self, item):
        return False

    def _get(self, items):
        if self.fnc_uses_lists:
            result = self._eval([items])
            return result[0]
        else:
            result = self._eval(items)
            return result

    def _get_list(self, items):
        if self.fnc_uses_lists:
            result = self._eval(items)
            return result
        else:
            return [self._eval(obj) for obj in items]

    def get_transformed_view(self, transform):
        def fnc(obj):
            return transform(self(obj))

        return fnc

class NODStore(NestableObjectDict):
    def __init__(self, name, key_class, dimensions, store, scope=None):
        super(NODStore, self).__init__()
        self.name = name
        self.dimensions = dimensions
        self.store = store
        self.key_class = key_class

        if scope is None:
            self.scope = self
        else:
            self.scope = scope

        self.max_save_buffer_size = 100

    def _add_new(self, items, values):
        [dict.__setitem__(self, item, value) for item, value in zip(items, values)]
        if len(self) > self.max_save_buffer_size:
            self.sync()

    @property
    def storage(self):
        return self.store.storage

    def sync(self):
        storable = { key : value for key, value in self.iteritems() if self.storage in key.idx }

        self.store.set_values(self.scope, storable.keys(), storable.values())
        # TODO: Allow to keep not saved values
        self.clear()

    def _get_key(self, item):
        if type(item) is tuple:
            if item[0].content_class is self.key_class:
                return item[1]
            else:
                return None

        elif self.storage in item.idx:
            return item.idx[self.storage]

        return None

    def _get(self, item):
        if item in self:
            return self[item]

        key = self._get_key(item)

        if key is None:
            return None

        return self._load(self, key)

    def _get_list(self, items):
        cached, missing = self._split_list_dict(self, items)

        keys = [self._get_key(item) for item in missing]
        replace = self._apply_some_list(self._load_list, keys)

        return self._replace_none(cached, replace)

    def _load(self, key):
        return self.store.get_value(self.scope, key)

    def _load_list(self, keys):
        return self.store.get_list_value(self.scope, keys)

    def _basetype(self, item):
        if type(item) is tuple:
            return item[0].content_class
        elif hasattr(item, 'base_cls'):
            return item.base_cls
        else:
            return type(item)

class NODMultiStore(NODStore):
    pass


class NODWrap(NestableObjectDict):
    def __init__(self, post):
        super(NODWrap, self).__init__()
        self.post = post

    def __getitem__(self, items):
        return self.post[items]

    def __setitem__(self, key, value):
        self.post[key] = value

class NODBuffer(NestableObjectDict):
    """
    Implements a dict with a buffer that reads sequentially ahead
    """

class NODCache(NestableObjectDict):
    """
    Implements a dict with intelligent caching of limited size
    """

@restores_as_stub_object
class OrderParameter(NODWrap):
    """
    Wrapper for a function that maps a snapshot to a number.

    Parameters
    ----------
    name : string
        A descriptive name of the orderparameter. It is used in the string
        representation.
    dimensions : int
        The number of dimensions of the output order parameter. So far this
        is not used and will be necessary or useful when storage is
        available
    storages : list of ConfigurationStorages()
        contains the list of storages that will be used.

    Attributes
    ----------
    name
    dimensions
    storages

    """
    def __init__(self, name, dimensions = 1):
        if type(name) is str and len(name) == 0:
            raise ValueError('name must be a non-empty string')

        self.pre_dict = NODTransform(self._pre_item)
        self.multi_dict = NODMultiStore()
        self.store_dict = NODStore(name, paths.Snapshot, dimensions, 'collectivevariable')
        self.cache_dict = NestableObjectDict()
        self.func_dict = NODFunction(None)
        self.func_dict._eval = self._eval

        self.name = name
        super(OrderParameter, self).__init__(
            post= self.func_dict + self.store_dict +
                  self.cache_dict + self.multi_dict + self.pre_dict
        )

    def _pre_item(self, items):
        item_type = self.store_dict._basetype(items)

        if item_type is  paths.Snapshot:
            return self._update(items)
        elif item_type is paths.Trajectory:
            return self._update(list(list.__iter__(items)))
        elif item_type is list:
            item_sub_type = self._basetype(items[0])
            if item_sub_type is paths.Snapshot:
                return self._update(items)
            else:
                return None
        else:
            return None

@restores_as_stub_object
class OP_RMSD_To_Lambda(OrderParameter):
    """
    Transforms the RMSD from `center` to a value between zero and one.

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
        trajectory = paths.Trajectory(items)
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
        trajectory = paths.Trajectory(items)

        # create an MDtraj trajectory out of it
        ptraj = trajectory.subset(self.atom_indices).md()

        # run the featurizer
        result = self.featurizer.partial_transform(ptraj)

        return result

@restores_as_stub_object
class OP_MD_Function(OrderParameter):
    """Make `OrderParameter` from `fcn` that takes mdtraj.trajectory as input.

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
        trajectory = paths.Trajectory(items)

        if self.topology is None:
            # first time ever compute the used topology for this orderparameter to construct the mdtraj objects
            self.topology = trajectory.topology.md

        t = trajectory.md(self.topology)
        return self.fcn(t, *args, **self.kwargs)

@restores_as_stub_object
class OP_Volume(OrderParameter):
    """ Make `Volume` into `OrderParameter`: maps to 0.0 or 1.0 """

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
    """Make any function `fcn` into an `OrderParameter`.

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

        trajectory = paths.Trajectory(items)

        return [ self._fcn(snap, *args, **self.kwargs) for snap in trajectory]
