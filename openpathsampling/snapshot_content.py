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

    def __init__(self, coordinates=None, box_vectors=None,
                 potential_energy=None):
        """
        Create a simulation configuration from either an OpenMM context or
        individually-specified components.

        Parameters
        ----------
        coordinates : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic coordinates (default: None)
        box_vectors : periodic box vectors (default: None)
            the periodic box vectors at current timestep (defautl: None)
        potential_energy : simtk.unit.Quantity of units energy/mole
            potential energy at current timestep (default: None)

        Attributes
        ----------
        coordinates : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic coordinates
        box_vectors : periodic box vectors
            the periodic box vectors
        potential_energy : simtk.unit.Quantity of units energy/mole
            potential energy
        idx : dict( Storage() : int )
            dict for storing the used index per storage
        topology : mdtraj.Topology()
            a reference to the used topology. This is necessary to allow
            export to mdtraj objects
        """

        super(Configuration, self).__init__()

        self.coordinates = None
        self.box_vectors = None
        self.potential_energy = None

        # TODO: Replace deepcopy by reference. Deepcopy is against immutable agreement
        if coordinates is not None:
            self.coordinates = copy.deepcopy(coordinates)
        if box_vectors is not None:
            self.box_vectors = copy.deepcopy(box_vectors)
        if potential_energy is not None:
            self.potential_energy = copy.deepcopy(potential_energy)

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
        '''
        Returns the number of atoms in the configuration
        '''
        return self.coordinates.shape[0]

    # =========================================================================
    # Utility functions
    # =========================================================================

    def copy(self, subset=None):
        """
        Returns a deep copy of the instance itself. If this object is saved
        it will not be stored as a separate object and consume additional
        memory. Should be avoided!

        Returns
        -------
        Configuration()
            the deep copy
        """

        if subset is None:
            this = Configuration(coordinates=self.coordinates,
                                 box_vectors=self.box_vectors,
                                 potential_energy=self.potential_energy
                                 )
        else:
            new_coordinates = self.coordinates[subset, :]
            # TODO: Keep old potential_energy? Is not correct but might be useful. Boxvectors are fine!
            this = Configuration(coordinates=new_coordinates,
                                 box_vectors=self.box_vectors,
                                 potential_energy=self.potential_energy
                                 )

        return this


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

    def __init__(self, velocities=None, kinetic_energy=None):
        """
        Create a simulation momentum from either an OpenMM context or
        individually-specified components.

        Parameters
        ----------
        velocities : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic velocities (default: None)
        kinetic_energy : simtk.unit.Quantity of units energy/mole
            kinetic energy at current timestep (default: None)

        Attributes
        ----------
        velocities : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic velocities
        kinetic_energy : simtk.unit.Quantity of units energy/mole
            kinetic energy
        idx : dict( Storage() : int )
            dict for storing the used index per storage
        """

        super(Momentum, self).__init__()

        self.velocities = None
        self.kinetic_energy = None

        if velocities is not None:
            self.velocities = copy.deepcopy(velocities)
        if kinetic_energy is not None:
            self.kinetic_energy = copy.deepcopy(kinetic_energy)

        return

    @property
    def n_atoms(self):
        '''
        Returns the number of atoms in the momentum
        '''
        return self.velocities.shape[0]

    # =========================================================================
    # Utility functions
    # =========================================================================

    def copy(self, subset=None, reversed=False):
        """
        Returns a deep copy of the instance itself. If this object will not
        be saved as a separate object and consumes additional memory. It is
        used to construct a reversed copy that can be stored or used to
        start a simulation. If the momentum is shallow it will be loaded for
        the copy

        Returns
        -------
        Momentum()
            the deep copy
        """

        if subset is None:
            new_velocities = self.velocities
        else:
            new_velocities = self.velocities[subset, :]
            # TODO: Keep old kinetic_energy? Is not correct but might be useful.

        if reversed:
            # Note the v *= -1.0 would be in place for numpy arrays. This here makes a copy!
            new_velocities = -1.0 * new_velocities

        this = Momentum(velocities=new_velocities, kinetic_energy=self.kinetic_energy)

        return this

    def reversed_copy(self, subset=None):
        """
        Create a copy and flips the velocities and erases the stored indices.
        If stores is will be treated as a new Momentum instance.
        Should be avoided.

        Returns
        -------
        Momentum()
            the deep copy with reversed velocities.
        """
        return self.copy(subset=subset, reversed=True)


