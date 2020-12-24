import copy

import numpy as np
from openpathsampling.netcdfplus import StorableObject, ObjectStore, WeakLRUCache
from openpathsampling.integration_tools import error_if_no_simtk_unit, unit

def unmask_quantity(quantity):
    """Force a maskedarray quantity to be unmasked.

    NetCDF keeps giving us masked arrays, even when we tell it not to.
    Masked arrays cause all kinds of havoc with other parts of the code.

    Parameters
    ----------
    quantity : simtk.unit.Quantity wrapping a numpy (masked) array
        quantity to unmask

    Returns
    -------
    simtk.unit.Quantity
        wraps a regular numpy array, not a masked array
    """
    try:
        q_unit = quantity.unit
    except AttributeError:
        # no units
        return quantity
    return np.array(quantity.value_in_unit(q_unit)) * q_unit


# =============================================================================
# SIMULATION CONFIGURATION
# =============================================================================

class StaticContainer(StorableObject):
    """
    Simulation configuration. Only Coordinates, the associated boxvectors
    and the potential_energy

    Parameters
    ----------
    coordinates : simtk.unit.Quantity wrapping Nx3 np array of dimension length
        atomic coordinates
    box_vectors : periodic box vectors
        the periodic box vectors
    engine : :class:`.DynamicsEngine`
        the engine that creating this data
    """

    # Class variables to store the global storage and the system context
    # describing the system to be saved as configuration_indices

    def __init__(self, coordinates, box_vectors, engine=None):
        super(StaticContainer, self).__init__()
        self.coordinates = copy.deepcopy(coordinates)
        self.box_vectors = copy.deepcopy(box_vectors)
        self.engine = engine

    @property
    def n_atoms(self):
        """
        Returns the number of atoms in the configuration
        """
        return self.coordinates.shape[0]

    def copy(self):
        """
        Returns a deep copy of the instance itself using a subset of coordinates.
        If this object is saved it will be stored as a separate object and
        consume additional memory.

        Returns
        -------
        Configuration()
            the reduced deep copy
        """

        return StaticContainer(coordinates=self.coordinates,
                               box_vectors=self.box_vectors,
                               engine=self.engine)

    def to_dict(self):
        # note: to_dict not used here in SimStore, so no need to change
        return {
            'coordinates': self.coordinates,
            'box_vectors': self.box_vectors
        }


class StaticContainerStore(ObjectStore):
    """
    An ObjectStore for Configuration. Allows to store Configuration() instances in a netcdf file.
    """
    def __init__(self):
        super(StaticContainerStore, self).__init__(StaticContainer, json=False)

    def to_dict(self):
        return {}

    def _save(self, configuration, idx):
        # Store configuration.
        self.vars['coordinates'][idx] = configuration.coordinates
        box_vectors = configuration.box_vectors

        if box_vectors is None:
            n_spatial = configuration.coordinates.shape[1]
            box_vectors = np.zeros((n_spatial, n_spatial))

        self.vars['box_vectors'][idx] = box_vectors

    def get(self, indices):
        return [self.load(idx) for idx in indices]

    def _load(self, idx):
        coordinates = self.vars["coordinates"][idx]
        box_vectors = self.vars["box_vectors"][idx]

        if not np.count_nonzero(box_vectors):
            box_vectors = None

        configuration = StaticContainer(coordinates=coordinates,
                                        box_vectors=box_vectors)

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

        variable = self.storage.variables[self.prefix + '_coordinates']

        return variable[frame_indices, atom_indices, :].astype(
            np.float32).copy()

    def initialize(self):
        super(StaticContainerStore, self).initialize()
        error_if_no_simtk_unit("StaticContainerStore")

        self.create_variable(
            'coordinates', 'numpy.float32',
            dimensions=('n_atoms', 'n_spatial'),
            description="coordinate of atom '{ix[1]}' in dimension " +
                        "'{ix[2]}' of configuration '{ix[0]}'.",
            chunksizes=('n_atoms', 'n_spatial'),
            simtk_unit=unit.nanometers)

        self.create_variable(
            'box_vectors', 'numpy.float32',
            dimensions=('n_spatial', 'n_spatial'),
            chunksizes=('n_spatial', 'n_spatial'),
            simtk_unit=unit.nanometers)


# =============================================================================
# SIMULATION MOMENTUM / VELOCITY
# =============================================================================

class KineticContainer(StorableObject):
    """
    Simulation momentum. Contains only velocities of all atoms and
    associated kinetic energies

    Attributes
    ----------
    velocities : simtk.unit.Quantity wrapping Nx3 np array of dimension length
        atomic velocities
    engine : :class:`.DynamicsEngine`
        the engine that creating this data
    """

    def __init__(self, velocities, engine=None):
        super(KineticContainer, self).__init__()
        self.velocities = copy.deepcopy(velocities)
        self.engine = engine

    def copy(self):
        """
        Returns a deep copy of the instance itself. If saved this object will
        be stored as a separate object and consume additional memory.

        Returns
        -------
        Momentum()
            the shallow copy
        """
        this = KineticContainer(velocities=self.velocities,
                                engine=self.engine)
        return this

    def to_dict(self):
        return {
            'velocities': self.velocities
        }


class KineticContainerStore(ObjectStore):
    """
    An ObjectStore for Momenta. Allows to store Momentum() instances in a netcdf file.
    """

    def __init__(self):
        super(KineticContainerStore, self).__init__(KineticContainer, json=False)

    def to_dict(self):
        return {}

    def _save(self, momentum, idx):
        self.vars['velocities'][idx, :, :] = momentum.velocities

    def _load(self, idx):
        velocities = self.vars['velocities'][idx]

        momentum = KineticContainer(velocities=velocities)
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

        v = self.variables['velocities']

        return v[frame_indices, atom_indices, :].astype(np.float32).copy()

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

    def initialize(self):
        """
        Initializes the associated storage to index momentums in it
        """
        error_if_no_simtk_unit("KineticContainerStore")

        super(KineticContainerStore, self).initialize()

        self.create_variable(
            'velocities', 'numpy.float32',
            dimensions=('n_atoms', 'n_spatial'),
            description="the velocity of atom 'atom' in dimension " +
                        "'coordinate' of momentum 'momentum'.",
            chunksizes=('n_atoms', 'n_spatial'),
            simtk_unit=unit.nanometers / unit.picoseconds)
