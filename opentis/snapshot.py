'''

@author: JD Chodera
@author: JH Prinz
'''

import copy

import numpy as np
import mdtraj as md


#=============================================================================================
# SIMULATION CONFIGURATION
#=============================================================================================
from wrapper import storable


@storable
class Configuration(object):

    """
    Simulation configuration. Only Coordinates, the associated boxvectors and the potential_energy

    """

    # Class variables to store the global storage and the system context describing the system to be safed as configuration_indices
    simulator = None
    load_lazy = True

    def __init__(self, context=None, simulator=None, coordinates=None, box_vectors=None, potential_energy=None, topology=None, idx=None):
        """
        Create a simulation configuration from either an OpenMM context or individually-specified components.

        Parameters
        ----------
        context : simtk.chem.openContext
            if not None, the current state will be queried to populate simulation configuration;
            otherwise, can specify individual components (default: None)
        simulator : Simulator()
            if not None, the context and the topology is taken from the simulator object. This
            should be the preferred way when using simulations
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
            a reference to the used topology. This is necessary to allow export to mdtraj objects

        """

        self._coordinates = None
        self._box_vectors = None
        self._potential_energy = None
        self.topology = None

        if simulator is not None:
            context = simulator.simulation.context
            self.topology = simulator.storage.topology

        if topology is not None:
            self.topology = topology

        if context is not None:
            # Get current state from OpenMM Context object.
            state = context.getState(getPositions=True, getEnergy=True)

            # Store the associated context
            self.context = context

            # Populate current configuration data.
            self._coordinates = state.getPositions(asNumpy=True)
            self._box_vectors = state.getPeriodicBoxVectors()
            self._potential_energy = state.getPotentialEnergy()
        else:
            if coordinates is not None: self._coordinates = copy.deepcopy(coordinates)
            if box_vectors is not None: self._box_vectors = copy.deepcopy(box_vectors)
            if potential_energy is not None: self._potential_energy = copy.deepcopy(potential_energy)

        if self._coordinates is not None:
            # Check for nans in coordinates, and raise an exception if something is wrong.
            if np.any(np.isnan(self._coordinates)):
                raise Exception("Some coordinates became 'nan'; simulation is unstable or buggy.")

        return

    @property
    def coordinates(self):
        if Configuration.load_lazy is True and self._coordinates is None and len(self.idx) > 0:
            # this uses the first storage and loads the velocities from there
            self.idx.iterkeys().next().configuration.update_coordinates(self)

        return self._coordinates

    @coordinates.setter
    def coordinates(self, value):
        if self._coordinates is None:
            self._coordinates = value
        else:
            raise ValueError("Cannot change coordinates once they are set")

    @property
    def box_vectors(self):
        if Configuration.load_lazy is True and self._box_vectors is None and len(self.idx) > 0:
            # this uses the first storage and loads the velocities from there
            self.idx.iterkeys().next().configuration.update_box_vectors(self)

        return self._box_vectors

    @box_vectors.setter
    def box_vectors(self, value):
        if self._box_vectors is None:
            self._box_vectors = value
        else:
            raise ValueError("Cannot change box_vector once they are set")

    @property
    def potential_energy(self):
        if Configuration.load_lazy is True and self._potential_energy is None and len(self.idx) > 0:
            # this uses the first storage and loads the velocities from there
            self.idx.iterkeys().next().configuration.update_potential_energy(self)

        return self._potential_energy

    @potential_energy.setter
    def potential_energy(self, value):
        if self._potential_energy is None:
            self._potential_energy = value
        else:
            raise ValueError("Cannot change potential_energy once they are set")

    def forget(self):
        """
        Will remove the stored coordinates from memory if they are stored in a file to save memory.
        Once the coordinates are accessed they are reloaded automatically
        """

        if Configuration.load_lazy and len(self.idx) > 0:
            self._coordinates = None
            self._box_vectors = None
            self._potential_energy = None

    #=============================================================================================
    # Comparison functions
    #=============================================================================================

    def __eq__(self, other):
        if self is other:
            return True
        for storage in self.idx:
            if storage in other.idx and other.idx[storage] == self.idx[storage]:
                return True

        return False

#        return self is other

#    def __hash__(self):
        # We need to make sure that a configuration from storage can be found. So just take the numpy
        # array of coordinates call a tostring and use this. That should be reasonably fast and should
        # only avoid checking all coordinates. The final check is really, if the elements were loaded
        # from the same idx in the same file (see __eq__)

#        return hash(self.coordinates.tostring())

    @property
    def atoms(self):
        '''
        Returns the number of atoms in the configuration
        '''
        return self.coordinates.shape[0]

    #=============================================================================================
    # Utility functions
    #=============================================================================================

    def copy(self):
        """
        Returns a deep copy of the instance itself. If this object is saved it will not be stored as a
        separate object and consume additional memory. Should be avoided!

        Returns
        -------
        Configuration()
            the deep copy
        """

        this = Configuration(coordinates=self.coordinates, box_vectors=self._box_vectors, potential_energy=self._potential_energy, topology=self.topology)
        return this

    def md(self):
        '''
        Returns a mdtraj.Trajectory() object that contains only one frame

        Returns
        -------
        mdtraj.Tractory
            the actual trajectory object. Can be used with all functions from mdtraj

        Notes
        -----
        Rather slow since the topology has to be made each time. Try to avoid it
        '''

        n_atoms = self.atoms

        output = np.zeros([1, n_atoms, 3], np.float32)
        output[0,:,:] = self.coordinates

        return md.Trajectory(output, self.topology)












#=============================================================================================
# SIMULATION MOMENTUM / VELOCITY
#=============================================================================================
@storable
class Momentum(object):
    """
    Simulation momentum. Contains only velocities of all atoms and associated kinetic energies
    """
    
    # Class variables to store the global storage and the system context describing the system to be safed as momentums
    simulator = None
    load_lazy = True

    def __init__(self, context=None, simulator=None, velocities=None, kinetic_energy=None, idx=None):
        """
        Create a simulation momentum from either an OpenMM context or individually-specified components.

        Parameters
        ----------
        context : simtk.chem.openContext
            if not None, the current state will be queried to populate simulation momentum; 
            otherwise, can specify individual components (default: None)
        simulator : Simulator()
            if not None, the context and the topology is taken from the simulator object. This
            should be the preferred way when using simulations
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
        
        self._velocities = None
        self._kinetic_energy = None

        if simulator is not None:
            context = simulator.simulation.context

        if context is not None:
            # Get current state from OpenMM Context object.
            state = context.getState(getVelocities=True, getEnergy=True)
            
            # Store the associated context
            self.context = context
            
            # Populate current momentum data.
            self._velocities = state.getVelocities(asNumpy=True)
            self._kinetic_energy = state.getKineticEnergy()
        else:
            if velocities is not None: self._velocities = copy.deepcopy(velocities)
            if kinetic_energy is not None: self._kinetic_energy = copy.deepcopy(kinetic_energy)

        # Check for nans in coordinates, and raise an exception if something is wrong.
#        if np.any(np.isnan(self.coordinates)):
#            raise Exception("Some coordinates became 'nan'; simulation is unstable or buggy.")

        return

    @property
    def velocities(self):
        if Momentum.load_lazy is True and self._velocities is None and len(self.idx) > 0:
            # this uses the first storage and loads the velocities from there
            self.idx.iterkeys().next().momentum.update_velocities(self)

        return self._velocities

    @velocities.setter
    def velocities(self, value):
        if self._velocities is None:
            self._velocities = value
        else:
            raise ValueError()

    @property
    def kinetic_energy(self):
        if Momentum.load_lazy is True and self._velocities is None and len(self.idx) > 0:
            # this uses the first storage and loads the velocities from there
            self.idx.iterkeys().next().momentum.update_kinetic_energy(self)

        return self._kinetic_energy

    @kinetic_energy.setter
    def kinetic_energy(self, value):
        if self._kinetic_energy is None:
            self._kinetic_energy = value
        else:
            raise ValueError()

    def forget(self):
        """
        Will remove the stored Momentum data from memory if they are stored in a file to save memory.
        Once the coordinates are accessed they are reloaded automatically
        """

        if Momentum.load_lazy and len(self.idx) > 0:
            self._velocities = None
            self._kinetic_energy = None

    @property
    def atoms(self):
        '''
        Returns the number of atoms in the momentum
        '''   
        return self.velocities.shape[0]

    #=============================================================================================
    # Utility functions
    #=============================================================================================

    def copy(self):
        """
        Returns a deep copy of the instance itself. If this object will not be saved as a
        separate object and consumes additional memory. It is used to construct a reversed
        copy that can be stored or used to start a simulation. If the momentum is shallow it
        will be loaded for the copy

        Returns
        -------
        Momentum()
            the deep copy
        """
        this = Momentum(velocities=self._velocities, kinetic_energy=self._kinetic_energy)
        this.idx = self.idx
        return this

    def reverse(self):
        """
        Flips the velocities and erases the stored indices. If stores is will be treated as a new Momentum instance.
        Should be avoided.

        """

        # This trick loads both, velocities and the kinetic energy. Otherwise we might run into trouble when removing the index
        self._velocities = -1.0 * self.velocities
        self.kinetic_energy
        self.idx = dict()

    
    def reversed_copy(self):
        """
        Create a copy and flips the velocities and erases the stored indices.
        If stores is will be treated as a new Momentum instance.
        Should be avoided.

        Returns
        -------
        Momentum()
            the deep copy with reversed velocities.
        """
        this = self.copy()
        this.reverse()
        return this












#=============================================================================================
# SIMULATION SNAPSHOT (COMPLETE FRAME WITH COORDINATES AND VELOCITIES)
#=============================================================================================

class Snapshot(object):
    """
    Simulation snapshot. Contains references to a configuration and momentum

    """
    
    # Class variables to store the global storage and the system context describing the system to be saved as snapshots
    # Hopefully these class member variables will not be needed any longer
    simulator = None

    def __init__(self, context=None, simulator=None, coordinates=None, velocities=None, box_vectors=None, potential_energy=None, kinetic_energy=None, configuration=None, momentum=None, reversed=False, topology=None):
        """
        Create a simulation snapshot from either an OpenMM context or individually-specified components.

        Parameters
        ----------
        context : simtk.chem.openContext
            if not None, the current state will be queried to populate simulation snapshot; 
            otherwise, can specify individual components (default: None)
        simulator : Simulator()
            if not None, the context and the topology is taken from the simulator object. This
            should be the preferred way when using simulations
        coordinates : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic coordinates (default: None)
        velocities : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic velocities (default: None)
        box_vectors : periodic box vectors (default: None)
            the periodic box vectors at current timestep (defautl: None)
        potential_energy : simtk.unit.Quantity of units energy/mole
            potential energy at current timestep (default: None)
        kinetic_energy : simtk.unit.Quantity of units energy/mole
            kinetic energy at current timestep (default: None)
            
        Attributes
        ----------
        coordinates : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic coordinates
        velocities : simtk.unit.Quantity wrapping Nx3 np array of dimension length
            atomic velocities
        box_vectors : periodic box vectors
            the periodic box vectors 
        potential_energy : simtk.unit.Quantity of units energy/mole
            potential energy
        kinetic_energy : simtk.unit.Quantity of units energy/mole
            kinetic energy
        idx : dict( Storage() : int )
            dict for storing the used index per storage
        
        """
        
        if configuration is None:
            self.configuration = Configuration()
        else:
            self.configuration = configuration
        if momentum is None:
            self.momentum = Momentum()
        else:
            self.momentum = momentum

        if simulator is not None:
            context = simulator.simulation.context
            self.configuration.topology = simulator.storage.topology

        if topology is not None:
            self.configuration.topology = topology

        self.reversed = reversed

        if context is not None:
            # Get current state from OpenMM Context object.
            state = context.getState(getPositions=True, getVelocities=True, getEnergy=True)
            
            # Store the associated context
            self.context = context
            
            # Populate current snapshot data.
            self.configuration._coordinates = state.getPositions(asNumpy=True)
            self.momentum._velocities = state.getVelocities(asNumpy=True)
            self.configuration._box_vectors = state.getPeriodicBoxVectors(asNumpy=True)
            self.configuration._potential_energy = state.getPotentialEnergy()
            self.momentum._kinetic_energy = state.getKineticEnergy()
        else:
            if coordinates is not None: self.configuration._coordinates = copy.deepcopy(coordinates)
            if velocities is not None: self.momentum._velocities = copy.deepcopy(velocities)
            if box_vectors is not None: self.configuration._box_vectors = copy.deepcopy(box_vectors)
            if potential_energy is not None: self.configuration._potential_energy = copy.deepcopy(potential_energy)
            if kinetic_energy is not None: self.momentum._kinetic_energy = copy.deepcopy(kinetic_energy)
            
        if self.configuration._coordinates is not None:
            # Check for nans in coordinates, and raise an exception if something is wrong.
            if np.any(np.isnan(self.configuration._coordinates)):
                raise Exception("Some coordinates became 'nan'; simulation is unstable or buggy.")
                
        pass

    @property
    def topology(self):
        """
        The mdtraj.Topology store in the configuration if present.
        """
        if self.configuration is not None:
            return self.configuration.topology
        return None
    
    @property
    def coordinates(self):
        """
        The coordinates in the configuration
        """
        return self.configuration.coordinates

    @property
    def velocities(self):
        """
        The velocities in the configuration. If the snapshot is reversed a copy of the original
        (unreversed) velocities is made which is then returned
        """
        if self.reversed:
            return self.momentum.reversed_copy().velocities
        else:
            return self.momentum.velocities
    
    @property
    def box_vectors(self):
        """
        The box_vectors in the configuration
        """
        return self.configuration.box_vectors
    
    @property
    def potential_energy(self):
        """
        The potential_energy in the configuration
        """
        return self.configuration.potential_energy
    
    @property
    def kinetic_energy(self):
        """
        The kinetic_energy in the momentum
        """
        return self.momentum.kinetic_energy
    
    @property
    def atoms(self):
        '''
        The number of atoms in the snapshot
        ''' 
        return self.coordinates.shape[0]

    @property
    def total_energy(self):
        '''
        The total energy (sum of potential and kinetic) of the snapshot
        '''   
        return self.kinetic_energy + self.potential_energy
    
    #=============================================================================================
    # Utility functions
    #=============================================================================================

    def copy(self):
        """
        Returns a shallow copy of the instance itself. The contained configuration and momenta are not
        copied.
        Returns
        -------
        Snapshot()
            the deep copy
        """
        this = Snapshot(configuration=self.configuration, momentum=self.momentum, reversed=self.reversed)
        return this
    
    def reversed_copy(self):
        """
        Returns a shallow reversed copy of the instance itself. The contained configuration and momenta are not
        copied and the momenta are marked reversed.
        Returns
        -------
        Snapshot()
            the deep copy
        """

        return self.copy().reverse()

    def reverse(self):
        """
        Reversed the momenta. This only flips a boolean and marks the given snapshot are reversed. This is fast and
        should be used instead of read velocity inversion.
        """

        self.reversed = not self.reversed
        return self
    
    def md(self):
        '''
        Returns a mdtraj Trajectory object that contains only one frame
        
        Notes
        -----        
        Rather slow since the topology has to be made each time. Try to avoid it
        '''        
        return self.configuration.md()
