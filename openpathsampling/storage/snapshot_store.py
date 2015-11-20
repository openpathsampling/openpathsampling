import numpy as np

from openpathsampling.snapshot import Snapshot, Configuration, Momentum
from openpathsampling.trajectory import Trajectory
from openpathsampling.netcdfplus import ObjectStore, LoaderProxy
from openpathsampling.tools import units_from_snapshot
from openpathsampling.netcdfplus import ExternalFile, StorableObject, lazy_loading_attributes

import mdtraj as md
import simtk.unit as u


class SnapshotStore(ObjectStore):
    """
    An ObjectStore for Snapshots in netCDF files.
    """

    def __init__(self):
        super(SnapshotStore, self).__init__(Snapshot, json=False)

    def _load(self, idx):
        """
        Load a snapshot from the storage.

        Parameters
        ----------
        idx : int
            the integer index of the snapshot to be loaded

        Returns
        -------
        snapshot : Snapshot
            the loaded snapshot instance
        """

        s_idx = int(idx / 2) * 2

        reversed_idx = 2 * s_idx + 1 - idx
        try:
            obj = self.cache[reversed_idx]
            snapshot = Snapshot(
                configuration=obj.configuration,
                momentum=obj.momentum,
                is_reversed=not obj.is_reversed,
                reversed_copy=LoaderProxy(self, reversed_idx)
            )
            return snapshot


        except KeyError:
            pass

        configuration = self.vars['configuration'][s_idx]
        momentum = self.vars['momentum'][s_idx]
        momentum_reversed = self.vars['momentum_reversed'][idx]

        snapshot = Snapshot(
            configuration=configuration,
            momentum=momentum,
            is_reversed=momentum_reversed,
            reversed_copy=LoaderProxy(self, reversed_idx)
        )

        return snapshot

    def _save(self, snapshot, idx):
        """
        Add the current state of the snapshot in the database.

        Parameters
        ----------
        snapshot : Snapshot()
            the snapshot to be saved
        idx : int or None
            if idx is not None the index will be used for saving in the storage.
            This might overwrite already existing trajectories!

        Notes
        -----
        This also saves all contained frames in the snapshot if not done yet.
        A single Snapshot object can only be saved once!
        """

        s_idx = int(idx / 2) * 2

        reversed_idx = 2 * s_idx + 1 - idx

        self.write('configuration', s_idx + 1, snapshot, to_lazy=False)
        self.vars['momentum'][s_idx + 1] = snapshot.momentum
        self.write('configuration', s_idx, snapshot)
        self.write('momentum', s_idx, snapshot)

        self.vars['momentum_reversed'][idx] = snapshot.is_reversed
        self.vars['momentum_reversed'][reversed_idx] = not snapshot.is_reversed

        reversed = snapshot._reversed
        snapshot._reversed = LoaderProxy(self, reversed_idx)
        reversed._reversed = LoaderProxy(self, idx)

        # mark reversed as stored
        self.index[reversed] = reversed_idx

    def _init(self):
        """
        Initializes the associated storage to index configuration_indices in it
        """
        super(SnapshotStore, self)._init()

        self.init_variable('configuration', 'lazyobj.configurations',
                           description="the snapshot index (0..n_configuration-1) of snapshot '{idx}'.",
                           chunksizes=(1,)
                           )

        self.init_variable('momentum', 'lazyobj.momenta',
                           description="the snapshot index (0..n_momentum-1) 'frame' of snapshot '{idx}'.",
                           chunksizes=(1,)
                           )

        self.init_variable('momentum_reversed', 'bool', chunksizes=(1,))

    # =============================================================================================
    # COLLECTIVE VARIABLE UTILITY FUNCTIONS
    # =============================================================================================

    @property
    def op_configuration_idx(self):
        """
        Returns aa function that returns for an object of this storage the idx

        Returns
        -------
        function
            the function that returns the idx of the configuration
        """

        def idx(obj):
            return self.index[obj.configuration]

        return idx

    @property
    def op_momentum_idx(self):
        """
        Returns aa function that returns for an object of this storage the idx

        Returns
        -------
        function
            the function that returns the idx of the configuration

        """

        def idx(obj):
            return self.index[obj.momentum]

        return idx

    def all(self):
        return Trajectory([LoaderProxy(self, idx) for idx in range(len(self))])


class MomentumStore(ObjectStore):
    """
    An ObjectStore for Momenta. Allows to store Momentum() instances in a netcdf file.
    """

    def __init__(self):
        super(MomentumStore, self).__init__(Momentum, json=False)

    def _save(self, momentum, idx):
        self.vars['velocities'][idx, :, :] = momentum.velocities

        if momentum.kinetic_energy is not None:
            self.vars['kinetic_energy'][idx] = momentum.kinetic_energy

    def _load(self, idx):
        velocities = self.vars['velocities'][idx]
        kinetic_energy = self.vars['kinetic_energy'][idx]

        momentum = Momentum(velocities=velocities, kinetic_energy=kinetic_energy)

        return momentum

    def velocities_as_numpy(self, frame_indices=None, atom_indices=None):
        """
        Return a block of stored velocities in the database as a numpy array.

        Parameters
        ----------
        frame_indices : list of int or None
            the indices of Momentum objects to be retrieved from the database.
            If `None` is specified then all indices are returned!
        atom_indices : list of int of None
            if not None only the specified atom_indices are returned. Might
            speed up reading a lot.
        """

        if frame_indices is None:
            frame_indices = slice(None)

        if atom_indices is None:
            atom_indices = slice(None)

        return self.variables['velocities'][frame_indices, atom_indices, :].astype(np.float32).copy()

    def velocities_as_array(self, frame_indices=None, atom_indices=None):
        """
        Returns a numpy array consisting of all velocities at the given indices

        Parameters
        ----------
        frame_indices : list of int
            momenta indices to be loaded
        atom_indices : list of int
            selects only the atoms to be returned. If None (Default) all atoms
            will be selected


        Returns
        -------
        numpy.ndarray, shape = (l,n)
            returns an array with `l` the number of frames and `n` the number
            of atoms
        """

        return self.velocities_as_numpy(frame_indices, atom_indices)

    def _init(self):
        """
        Initializes the associated storage to index momentums in it
        """

        super(MomentumStore, self)._init()

        n_atoms = self.storage.n_atoms
        n_spatial = self.storage.n_spatial

        units = units_from_snapshot(self.storage._template)

        self.init_variable('velocities', 'numpy.float32',
                           dimensions=('atom', 'spatial'),
                           description="the velocity of atom 'atom' in dimension " +
                                       "'coordinate' of momentum 'momentum'.",
                           chunksizes=(1, n_atoms, n_spatial),
                           simtk_unit=units['velocity']
                           )

        self.init_variable('kinetic_energy', 'float',
                           chunksizes=(1,),
                           simtk_unit=units['energy']
                           )


class ConfigurationStore(ObjectStore):
    def __init__(self):
        super(ConfigurationStore, self).__init__(Configuration, json=False)

    def _save(self, configuration, idx):
        # Store configuration.
        self.vars['coordinates'][idx] = configuration.coordinates

        if configuration.potential_energy is not None:
            self.vars['potential_energy'][idx] = configuration.potential_energy

        if configuration.box_vectors is not None:
            self.vars['box_vectors'][idx] = configuration.box_vectors

    def get(self, indices):
        return [self.load(idx) for idx in indices]

    def _load(self, idx):
        coordinates = self.vars["coordinates"][idx]
        box_vectors = self.vars["box_vectors"][idx]
        potential_energy = self.vars["potential_energy"][idx]

        configuration = Configuration(coordinates=coordinates, box_vectors=box_vectors,
                                      potential_energy=potential_energy)
        configuration.topology = self.storage.topology

        return configuration

    def coordinates_as_numpy(self, frame_indices=None, atom_indices=None):
        """
        Return the atom coordinates in the storage for given frame indices
        and atoms

        Parameters
        ----------
        frame_indices : list of int or None
            the frame indices to be included. If None all frames are returned
        atom_indices : list of int or None
            the atom indices to be included. If None all atoms are returned

        Returns
        -------
        numpy.array, shape=(n_frames, n_atoms)
            the array of atom coordinates in a float32 numpy array

        """
        if frame_indices is None:
            frame_indices = slice(None)

        if atom_indices is None:
            atom_indices = slice(None)

        return self.storage.variables[self.prefix + '_coordinates'][frame_indices, atom_indices, :].astype(
            np.float32).copy()

    def _init(self):
        super(ConfigurationStore, self)._init()
        n_atoms = self.storage.n_atoms
        n_spatial = self.storage.n_spatial

        units = units_from_snapshot(self.storage._template)

        self.init_variable('coordinates', 'numpy.float32',
                           dimensions=('atom', 'spatial'),
                           description="coordinate of atom '{ix[1]}' in dimension " +
                                       "'{ix[2]}' of configuration '{ix[0]}'.",
                           chunksizes=(1, n_atoms, n_spatial),
                           simtk_unit=units['length']
                           )

        self.init_variable('box_vectors', 'numpy.float32',
                           dimensions=('spatial', 'spatial'),
                           chunksizes=(1, n_spatial, n_spatial),
                           simtk_unit=units['length']
                           )

        self.init_variable('potential_energy', 'float',
                           chunksizes=(1,),
                           simtk_unit=units['energy']
                           )


class ExternalMDConfigurationStore(ObjectStore):

    def __init__(self):
        super(ExternalMDConfigurationStore, self).__init__(ExternalMDConfiguration, json=False)

    def _save(self, mdframe, idx):
        self.write('mdfile', idx, mdframe)
        self.vars['frame'][idx] = mdframe.frame

    def _load(self, idx):
        mdfile = self.vars['mdfile'][idx]
        frame = self.vars['frame'][idx]

        return ExternalMDConfiguration(mdfile, frame, self.storage.topology)

    def _init(self):
        super(ExternalMDConfigurationStore, self)._init()

        self.init_variable('mdfile', 'lazyobj.mdfiles',
                           chunksizes=(1, ),
                           )

        self.init_variable('frame', 'int',
                           chunksizes=(1, ),
                           )

class ExternalMDFile(ExternalFile):
    """A Reference to an external file
    """

    def __init__(self, url):
        super(ExternalMDFile, self).__init__(url)
        self._mdfile = None

    @property
    def file(self):
        if self._mdfile is None:
            self._mdfile = md.load(self.url)

        return self._mdfile

@lazy_loading_attributes('mdfile')
class ExternalMDConfiguration(StorableObject):
    """A Reference to an external file
    """

    def __init__(self, mdfile, frame, topology=None):
        super(ExternalMDConfiguration, self).__init__()
        self.mdfile = mdfile
        self.frame = frame
        self.topology = topology

    @property
    def configuration(self):
        mdtrajectory = self.mdfile.file

        coord = u.Quantity(mdtrajectory.xyz[self.frame], u.nanometers)
        if mdtrajectory.unitcell_vectors is not None:
            box_v = u.Quantity(mdtrajectory.unitcell_vectors[self.frame],
                               u.nanometers)
        else:
            box_v = None

        configuration = Configuration(
            coordinates=coord,
            box_vectors=box_v,
            potential_energy=u.Quantity(0.0, u.kilojoule_per_mole),
            topology=self.topology
        )

        return configuration

class MultiConfigurationStore(ObjectStore):
    """
    A store that represents Configurations from multiple sources
    """

    def __init__(self):
        super(MultiConfigurationStore, self).__init__(Configuration, json=False)

    def _load(self, idx):
        # try loading from internal first
        configuration = self.vars['intconf'][idx]

        if configuration is None:
            # try external
            configuration = self.vars['extconf'][idx].configuration

        return configuration

    def persist(self, idx):
        int_idx = self.variables['intconf'][idx]
        ext_idx = self.variables['extconf'][idx]
        if int_idx == -1 and ext_idx >= 0:
            self.vars['intconf'][idx] = self.vars['extconf'][idx].configuration

    def _save(self, obj, idx):
        # saving is only done internal
        # this function is never called by proxies

        self.vars['intconf'][idx] = obj
        self.vars['extconf'][idx] = None

    def save_as_external(self, mdfile, frame):
        obj = ExternalMDConfiguration(mdfile, frame)
        idx = len(self)
        self.vars['intconf'][idx] = None
        self.vars['extconf'][idx] = obj

        return idx

    def _init(self):
        """
        Initializes the associated storage to index configuration_indices in it
        """
        super(MultiConfigurationStore, self)._init()

        self.init_variable('intconf', 'obj.intconf',
                           description="the snapshot index (0..n_configuration-1) of snapshot '{idx}'.",
                           chunksizes=(1,)
                           )

        self.init_variable('extconf', 'obj.extconf',
                           chunksizes=(1,)
                           )
