'''

@author: JD Chodera
@author: JH Prinz
'''

import copy
import numpy as np
import mdtraj as md

from simtk.unit import nanosecond, picosecond, nanometers, nanometer, picoseconds, femtoseconds, femtosecond, kilojoules_per_mole, Quantity

#=============================================================================================
# SIMULATION SNAPSHOT 
#=============================================================================================

class Configuration(object):

    """
    Simulation configuration. Only Coordinates, the associated boxvectors and the potential_energy

    """
    
    # Class variables to store the global storage and the system context describing the system to be safed as configuration_indices
    storage = None
    simulator = None
    
    def __init__(self, context=None, coordinates=None, box_vectors=None, potential_energy=None, simulator=None, topology=None):
        """
        Create a simulation configuration from either an OpenMM context or individually-specified components.

        Parameters
        ----------
        context : simtk.chem.openContext
            if not None, the current state will be queried to populate simulation configuration; 
            otherwise, can specify individual components (default: None)
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
        idx : int
            index for storing in the storage or for using with caching
        
        """
        
        self.idx = dict()        # potential idx in a netcdf storage, if 0 then not stored yet. Attention! Cannot be stored in 2 repositories at the same time
        self.coordinates = None
        self.box_vectors = None
        self.potential_energy = None
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
            self.coordinates = state.getPositions(asNumpy=True)
            self.box_vectors = state.getPeriodicBoxVectors()
            self.potential_energy = state.getPotentialEnergy()
        else:
            if coordinates is not None: self.coordinates = copy.deepcopy(coordinates)
            if box_vectors is not None: self.box_vectors = copy.deepcopy(box_vectors)
            if potential_energy is not None: self.potential_energy = copy.deepcopy(potential_energy)

        if self.coordinates is not None:
            # Check for nans in coordinates, and raise an exception if something is wrong.
            if np.any(np.isnan(self.coordinates)):
                raise Exception("Some coordinates became 'nan'; simulation is unstable or buggy.")

        return

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

    def __hash__(self):
        # We need to make sure that a configuration from storage can be found. So just take the numpy
        # array of coordinates call a tostring and use this. That should be reasonably fast and should
        # only avoid checking all coordinates. The final check is really, if the elements were loaded
        # from the same idx in the same file (see __eq__)
        return hash(self.coordinates.tostring())

    @property
    def atoms(self):
        '''
        The number of atoms in the configuration
        '''   
        return self.coordinates.shape[0]
    
    #=============================================================================================
    # Utility functions
    #=============================================================================================

    def copy(self):
        this = Configuration(coordinates=self.coordinates, box_vectors=self.box_vectors, potential_energy=self.potential_energy, topology=self.topology)
        return this

    def md(self):
        '''
        Returns a mdtraj Trajectory object that contains only one frame
        
        Notes
        -----        
        Rather slow since the topology has to be made each time. Try to avoid it
        '''        
        
        n_atoms = self.atoms
                            
        output = np.zeros([1, n_atoms, 3], np.float32)
        output[0,:,:] = self.coordinates

        return md.Trajectory(output, self.topology)


class Momentum(object):
    """
    Simulation momentum. Contains only velocities of all atoms and associated kinetic energies
    """
    
    # Class variables to store the global storage and the system context describing the system to be safed as momentums
    storage = None
    simulator = None
    
    def __init__(self, context=None, coordinates=None, velocities=None, box_vectors=None, potential_energy=None, kinetic_energy=None):
        """
        Create a simulation momentum from either an OpenMM context or individually-specified components.

        Parameters
        ----------
        context : simtk.chem.openContext
            if not None, the current state will be queried to populate simulation momentum; 
            otherwise, can specify individual components (default: None)
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
        idx : int
            index for storing in the storage or for using with caching
        
        """
        
        self.idx = dict()        # potential idx in a netcdf storage, if 0 then not stored yet. Attention! Cannot be stored in 2 repositories at the same time
        self.velocities = None
        self.kinetic_energy = None
        
        if context is not None:
            # Get current state from OpenMM Context object.
            state = context.getState(getVelocities=True, getEnergy=True)
            
            # Store the associated context
            self.context = context
            
            # Populate current momentum data.
            self.velocities = state.getVelocities(asNumpy=True)
            self.kinetic_energy = state.getKineticEnergy()
        else:
            if velocities is not None: self.velocities = copy.deepcopy(velocities)
            if kinetic_energy is not None: self.kinetic_energy = copy.deepcopy(kinetic_energy)                       

        # Check for nans in coordinates, and raise an exception if something is wrong.
#        if np.any(np.isnan(self.coordinates)):
#            raise Exception("Some coordinates became 'nan'; simulation is unstable or buggy.")

        return
    
    @property
    def atoms(self):
        '''
        The number of atoms in the momentum
        '''   
        return self.velocities.shape[0]

    #=============================================================================================
    # Utility functions
    #=============================================================================================

    def copy(self):
        this = Momentum(velocities=self.velocities, kinetic_energy=self.kinetic_energy)
        return this

    def reverse(self):
        self.velocities *= -1.0
    
    def reversed_copy(self):
        this = self.copy()
        this.reverse()
        return this
     
    #=============================================================================================
    # Storage functions
    #=============================================================================================
    

class Snapshot(object):
    """
    Simulation snapshot.

    """
    
    # Class variables to store the global storage and the system context describing the system to be safed as snapshots
    storage = None
    simulator = None
    
    def __init__(self, context=None, simulator=None, coordinates=None, velocities=None, box_vectors=None, potential_energy=None, kinetic_energy=None, configuration=None, momentum=None, reversed=False, topology=None):
        """
        Create a simulation snapshot from either an OpenMM context or individually-specified components.

        Parameters
        ----------
        context : simtk.chem.openContext
            if not None, the current state will be queried to populate simulation snapshot; 
            otherwise, can specify individual components (default: None)
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
        idx : int
            index for storing in the storage or for using with caching
        
        """
        
#        self.idx = 0        # potential idx in a netcdf storage, if 0 then not stored yet. Attention! Cannot be stored in 2 repositories at the same time
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
            self.configuration.coordinates = state.getPositions(asNumpy=True)
            self.momentum.velocities = state.getVelocities(asNumpy=True)
            self.configuration.box_vectors = state.getPeriodicBoxVectors(asNumpy=True)
            self.configuration.potential_energy = state.getPotentialEnergy()
            self.momentum.kinetic_energy = state.getKineticEnergy()            
        else:
            if coordinates is not None: self.configuration.coordinates = copy.deepcopy(coordinates)
            if velocities is not None: self.momentum.velocities = copy.deepcopy(velocities)
            if box_vectors is not None: self.configuration.box_vectors = copy.deepcopy(box_vectors)
            if potential_energy is not None: self.configuration.potential_energy = copy.deepcopy(potential_energy)
            if kinetic_energy is not None: self.momentum.kinetic_energy = copy.deepcopy(kinetic_energy)                       
            
        if self.coordinates is not None:
            # Check for nans in coordinates, and raise an exception if something is wrong.
            if np.any(np.isnan(self.coordinates)):
                raise Exception("Some coordinates became 'nan'; simulation is unstable or buggy.")
                
        pass

    @property
    def topology(self):
        if self.configuration is not None:
            return self.configuration.topology
        return None
    
    @property
    def coordinates(self):
        return self.configuration.coordinates

    @property
    def velocities(self):
        if self.reversed:
            return self.momentum.reversed_copy().velocities
        else:
            return self.momentum.velocities
    
    @property
    def box_vectors(self):
        return self.configuration.box_vectors
    
    @property
    def potential_energy(self):
        return self.configuration.potential_energy
    
    @property
    def kinetic_energy(self):
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
        this = Snapshot(configuration=self.configuration, momentum=self.momentum, reversed=self.reversed)
        return this
    
    def reversed_copy(self):
        return self.copy().reverse()

    def reverse(self):
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
    
    def md_topology(self): 
        '''
        Returns a mdtraj topology object that can be used with the stored snapshot
        '''   
        return self.configuration.md_topology()
    #=============================================================================================
    # Storage functions
    #=============================================================================================
