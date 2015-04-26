###############################################################
# | CLASS Order Parameter
###############################################################

import mdtraj as md
import openpathsampling as paths
import chaindict as cd
import collections
from openpathsampling.todict import ops_object


@ops_object
class CollectiveVariable(cd.Wrap):
    """
    Wrapper for a function that maps a snapshot to a number.

    Parameters
    ----------
    name : string
        A descriptive name of the collectivevariable. It is used in the string
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

    def __init__(self, name, dimensions=1):
        if (type(name) is not str and type(name) is not unicode) or len(name) == 0:
            print type(name), len(name)
            raise ValueError('name must be a non-empty string')

        self.name = name

        self.single_dict = cd.ExpandSingle()
        self.pre_dict = cd.Transform(self._pre_item)
        self.multi_dict = cd.ExpandMulti()
        self.store_dict = cd.MultiStore('collectivevariables', name,
                                            dimensions, self)
        self.cache_dict = cd.ChainDict()
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
                key : value for key, value in self.cache_dict
                    if type(key) is tuple or storage in key.idx
            }
            self.store_dict.cod_stores[storage].post.update(stored)
            self.store_dict.cod_stores[storage].update(stored)

    def sync(self, storage):
        """
        Sync this collectivevariable with attached storages

        Parameters
        ----------
        store : CollectiveVariableStore or None
            the store to be used, otherwise all underlying storages are synced
        store_cache : bool
            if `False` (default) the store will only store new values. If
            `True` also the cached values will be stored. This is much more
            costly and is usually only run once, when the collectivevariable is
            saved the first time
        """
        self.store_dict.update_nod_stores()
        if storage in self.store_dict.cod_stores:
            self.store_dict.cod_stores[storage].sync()

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

    _compare_keys = ['name']

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

@ops_object
class CV_Featurizer(CollectiveVariable):
    """
    An CollectiveVariable that uses an MSMBuilder3 featurizer as the logic

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
        super(CV_Featurizer, self).__init__(name,
                                            dimensions=featurizer.n_features)

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


@ops_object
class CV_MD_Function(CollectiveVariable):
    """Make `CollectiveVariable` from `fcn` that takes mdtraj.trajectory as input.

    Examples
    -------
    >>> # To create an order parameter which calculates the dihedral formed
    >>> # by atoms [7,9,15,17] (psi in Ala dipeptide):
    >>> import mdtraj as md
    >>> psi_atoms = [7,9,15,17]
    >>> psi_orderparam = CV_Function("psi", md.compute_dihedrals,
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
        super(CV_MD_Function, self).__init__(name)
        self.fcn = fcn
        self.kwargs = kwargs
        self.topology = None
        return

    def _eval(self, items, *args):
        trajectory = paths.Trajectory(items)

        if self.topology is None:
            # first time ever compute the used topology for this collectivevariable to construct the mdtraj objects
            self.topology = trajectory.topology.md

        t = trajectory.md(self.topology)
        return self.fcn(t, *args, **self.kwargs)


@ops_object
class CV_Volume(CollectiveVariable):
    """ Make `Volume` into `CollectiveVariable`: maps to 0.0 or 1.0 """

    def __init__(self, name, volume):
        """
        Parameters
        ----------

        """

        super(CV_Volume, self).__init__(name)
        self.volume = volume

    def _eval(self, items):
        result = [float(self.volume(item)) for item in items]
        return result

    def to_dict(self):
        return {
            'name' : self.name,
            'volume' : self.volume,
        }

    @staticmethod
    def from_dict(self, dct):
        return CV_Function(
            name=dct['name'],
            volume=dct['volume']
        )

    def __eq__(self, other):
        """Override the default Equals behavior"""
        if isinstance(other, self.__class__):
            if self.name != other.name:
                return False
            if self._fcn.func_code.op_code != other._fcn.func_code.op_code:
                # Compare Bytecode. Not perfect, but should be good enough
                return False

            return True

        return NotImplemented

    def __ne__(self, other):
        """Define a non-equality test"""
        if isinstance(other, self.__class__):
            return not self.__eq__(other)
        return NotImplemented

@ops_object
class CV_Function(CollectiveVariable):
    """Make any function `fcn` into an `CollectiveVariable`.

    Examples
    -------
    >>> # To create an order parameter which calculates the dihedral formed
    >>> # by atoms [7,9,15,17] (psi in Ala dipeptide):
    >>> import mdtraj as md
    >>> psi_atoms = [6,8,14,16]
    >>> psi_orderparam = CV_Function("psi", md.compute_dihedrals,
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
        super(CV_Function, self).__init__(name)
        self._fcn = fcn
        self.kwargs = kwargs
        return

    def to_dict(self):
        return {
            'name' : self.name,
            'fcn' : self.fcn,
            'kwargs' : self.kwargs
        }

    @staticmethod
    def from_dict(self, dct):
        return CV_Function(
            name=dct['name'],
            fcn=dct['fcn'],
            **dct['kwargs']
        )

    def __eq__(self, other):
        """Override the default Equals behavior"""
        if isinstance(other, self.__class__):
            if self.name != other.name:
                return False
            if self._fcn.func_code.op_code != other._fcn.func_code.op_code:
                # Compare Bytecode. Not perfect, but should be good enough
                return False

            return True

        return NotImplemented

    def __ne__(self, other):
        """Define a non-equality test"""
        if isinstance(other, self.__class__):
            return not self.__eq__(other)
        return NotImplemented

    def _eval(self, items, *args):
        trajectory = paths.Trajectory(items)

        return [self._fcn(snap, *args, **self.kwargs) for snap in trajectory]
