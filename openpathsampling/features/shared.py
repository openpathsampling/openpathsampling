import copy

import numpy as np
from simtk import unit as u

from openpathsampling.netcdfplus import StorableObject


# =============================================================================
# SIMULATION CONFIGURATION
# =============================================================================

class Configuration(StorableObject):
    """
    Simulation configuration. Only Coordinates, the associated boxvectors
    and the potential_energy
    """

    # Class variables to store the global storage and the system context
    # describing the system to be safed as configuration_indices

    def __init__(self, coordinates=None, box_vectors=None):
        """
        Create a simulation configuration from either an OpenMM context or
        individually-specified components.

        Parameters
        ----------
        coordinates : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic coordinates (default: None)
        box_vectors : periodic box vectors (default: None)
            the periodic box vectors at current timestep (defautl: None)

        Attributes
        ----------
        coordinates : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic coordinates
        box_vectors : periodic box vectors
            the periodic box vectors
        """

        super(Configuration, self).__init__()

        self.coordinates = None
        self.box_vectors = None

        # TODO: Replace deepcopy by reference. Deepcopy is against immutable agreement
        if coordinates is not None:
            self.coordinates = copy.deepcopy(coordinates)
        if box_vectors is not None:
            self.box_vectors = copy.deepcopy(box_vectors)

        if self.coordinates is not None:
            # Check for nans in coordinates, and raise an exception if
            # something is wrong.
            if type(self.coordinates) is u.Quantity:
                coords = self.coordinates._value
            else:
                coords = self.coordinates

            if np.any(np.isnan(coords)):
                raise ValueError(
                    "Some coordinates became 'nan'; simulation is unstable or buggy.")

        return

    # =========================================================================
    # Comparison functions
    # =========================================================================

    def __eq__(self, other):
        if self is other:
            return True

            # This is not good since this code requires knowledge about storage
            # I remove it since it is not used yet anyway
            # If we want to figure out if two Snapshots are loaded from two different
            # instances of the storage we should put this logic into storages

        #        for storage in self.idx:
        #            if storage in other.idx and other.idx[storage] == self.idx[storage]:
        #                return True

        return False

    @property
    def n_atoms(self):
        """
        Returns the number of atoms in the configuration
        """
        return self.coordinates.shape[0]

    # =========================================================================
    # Utility functions
    # =========================================================================

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

        # TODO: Keep old potential_energy? Is not correct but might be useful. Boxvectors are fine!
        return Configuration(coordinates=self.coordinates,
                             box_vectors=self.box_vectors
                             )


# =============================================================================
# SIMULATION MOMENTUM / VELOCITY
# =============================================================================

class Momentum(StorableObject):
    """
    Simulation momentum. Contains only velocities of all atoms and
    associated kinetic energies
    """

    # Class variables to store the global storage and the system context
    # describing the system to be safed as momentums

    def __init__(self, velocities=None):
        """
        Create a simulation momentum from either an OpenMM context or
        individually-specified components.

        Parameters
        ----------
        velocities : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic velocities (default: None)

        Attributes
        ----------
        velocities : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic velocities
        """

        super(Momentum, self).__init__()

        self.velocities = None

        if velocities is not None:
            self.velocities = copy.deepcopy(velocities)

        return

    @property
    def n_atoms(self):
        """
        Returns the number of atoms in the momentum
        """
        return self.velocities.shape[0]

    # =========================================================================
    # Utility functions
    # =========================================================================

    def copy(self):
        """
        Returns a deep copy of the instance itself. If saved this object will
        be stored as a separate object and consume additional memory.

        Returns
        -------
        Momentum()
            the shallow copy
        """

        this = Momentum(velocities=self.velocities)

        return this
