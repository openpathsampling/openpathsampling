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
    
    # Class variables to store the global storage and the system context describing the system to be safed as configurations
    storage = None
    simulator = None
    
    def __init__(self, context=None, coordinates=None, box_vectors=None, potential_energy=None):
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
        
        self.idx = 0        # potential idx in a netcdf storage, if 0 then not stored yet. Attention! Cannot be stored in 2 repositories at the same time
        self.coordinates = None
        self.box_vectors = None
        self.potential_energy = None
        
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

    @property
    def atoms(self):
        '''
        The number of atoms in the configuration
        '''   
        return len(self.coordinates.shape[0])  
    
    #=============================================================================================
    # Utility functions
    #=============================================================================================
    
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
        
        topology = self.md_topology()
                                                         
        return md.Trajectory(output, topology)      

    
    def md_topology(self): 
        '''
        Returns a mdtraj topology object that can be used with the stored configuration
        '''   
        return md.Topology.from_openmm(Configuration.simulator.simulation.topology)
 
    #=============================================================================================
    # Storage functions
    #=============================================================================================
    
    def save(self, idx = None):
        """
        Save positions, velocities, boxvectors and energies of current iteration to NetCDF file.
        
        Notes
        -----        
        We need to allow for reversed configurations to save memory. Would be nice        
        """
        
        storage = Configuration.storage
        
        if self.idx == 0:
            if idx is None:
                idx = Configuration.load_free()

            # Store configuration.
            storage.ncfile.variables['configuration_coordinates'][idx,:,:] = (self.coordinates / nanometers).astype(np.float32)
            if self.potential_energy is not None: 
                storage.ncfile.variables['configuration_potential'][idx] = self.potential_energy / kilojoules_per_mole                                
#            storage.ncfile.variables['configuration_box_vectors'][idx,:] = (self.box_vectors / nanometers).astype(np.float32)
            
            # store ID# for later reference in configuration object
            self.idx = idx
    
            # Force sync to disk to avoid data loss.
            storage.ncfile.sync()

        return

    @staticmethod
    def get(indices):
        return [Configuration.load(idx) for idx in indices ]
        

    @staticmethod
    def load_number():
        '''
        Load the number of stored configurations
        
        Returns
        -------
        number (int) - number of stored configurations
        '''
        length = int(len(Configuration.storage.ncfile.dimensions['configuration'])) - 1
        if length < 0:
            length = 0
        return length
    
    @staticmethod
    def load_free():
        '''
        Return the number of the next free ID
        '''
        return Configuration.load_number() + 1
    
    @staticmethod
    def load(idx):
        '''
        Load a configuration from the storage
        
        Parameters
        ----------
        idx : int
            index of the configuration in the database 'idx' > 0
        
        Returns
        -------
        configuration : configuration
            the configuration
        '''
        
        storage = Configuration.storage
        
        #TODO: Check, for some reason some idx are given as numpy.in32 and netcdf4 is not compatible with indices given in this format!!!!!
        idx = int(idx)
        
        x = storage.ncfile.variables['configuration_coordinates'][idx,:,:].astype(np.float32).copy()
        coordinates = Quantity(x, nanometers)                
        b = storage.ncfile.variables['configuration_box_vectors'][idx]
        box_vectors = Quantity(b, nanometers)              
        V = storage.ncfile.variables['configuration_potential'][idx]
        potential_energy = Quantity(V, kilojoules_per_mole)
    
        configuration = Configuration(coordinates=coordinates, box_vectors = box_vectors, potential_energy=potential_energy)
        configuration.idx = idx

        return configuration

    @staticmethod
    def coordinates_as_numpy(self, frame_indices=None, atom_indices=None):
        
        if frame_indices is None:
            frame_indices = slice(None)
        
        if atom_indices is None:
            atom_indices = slice(None)

        return self.ncfile.variables['configuration_coordinates'][frame_indices,atom_indices,:].astype(np.float32).copy()
    
    @staticmethod
    def _restore_netcdf(storage):
        """
        Fill in missing part after the storage has been loaded from a file and is not initialize freshly
        
        """
                
        Configuration.storage = storage
                
        return
        
    @staticmethod
    def _init_netcdf(storage):
        '''
        Initializes the associated storage to save configurations in it
        '''           
        # save associated storage in class variable for all configuration instances to access
        
#        ncgrp = storage.ncfile.createGroup('configuration')
        
        Configuration.storage = storage
        ncgrp = storage.ncfile
        
        system = Snapshot.simulator.simulation.system
        
        # define dimensions used in configurations
        ncgrp.createDimension('configuration', 0)                       # unlimited number of configurations
        if 'atom' not in ncgrp.dimensions:
            ncgrp.createDimension('atom', system.getNumParticles())    # number of atoms in the simulated system
        
        if 'spatial' not in ncgrp.dimensions:
            ncgrp.createDimension('spatial', 3) # number of spatial dimensions        

        # define variables for configurations
        ncvar_configuration_coordinates          = ncgrp.createVariable('configuration_coordinates', 'f', ('configuration','atom','spatial'))
        ncvar_configuration_box_vectors          = ncgrp.createVariable('configuration_box_vectors', 'f', ('configuration', 'spatial'))
        ncvar_configuration_potential            = ncgrp.createVariable('configuration_potential',   'f', ('configuration'))

        # Define units for configuration variables.
        setattr(ncvar_configuration_coordinates, 'units', 'nm')
        setattr(ncvar_configuration_box_vectors,  'units', 'nm')
        setattr(ncvar_configuration_potential,   'units', 'kJ/mol')
        
        # Define long (human-readable) names for variables.
        setattr(ncvar_configuration_coordinates,   "long_name", "coordinates[configuration][atom][coordinate] are coordinate of atom 'atom' in dimension 'coordinate' of configuration 'configuration'.")

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
        
        self.idx = 0        # potential idx in a netcdf storage, if 0 then not stored yet. Attention! Cannot be stored in 2 repositories at the same time
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
        return len(self.coordinates.shape[0])

    #=============================================================================================
    # Utility functions
    #=============================================================================================

    def reverse(self):
        self.velocities *= 1.0
    
    def reversed_copy(self):
        this = copy.deepcopy(self)
        this.reverse()
        return this
     
    #=============================================================================================
    # Storage functions
    #=============================================================================================
    
    def save(self, idx = None):
        """
        Save velocities and kinetic energies of current iteration to NetCDF file.
        """
        
        storage = Momentum.storage
        
        if self.idx == 0:
            if idx is None:
                idx = Momentum.load_free()

            # Store momentum.
            storage.ncfile.variables['momentumvelocities'][idx,:,:] = (self.velocities / (nanometers / picoseconds)).astype(np.float32)
            if self.kinetic_energy is not None:
                storage.ncfile.variables['momentum_kinetic'][idx] = self.kinetic_energy / kilojoules_per_mole
            
            # store ID# for later reference in Momentum object
            self.idx = idx
    
            # Force sync to disk to avoid data loss.
            storage.ncfile.sync()

        return

    @staticmethod
    def load_number():
        '''
        Load the number of stored momentums
        
        Returns
        -------
        number (int) - number of stored momentums
        '''
        length = int(len(Momentum.storage.ncfile.dimensions['momentum'])) - 1
        if length < 0:
            length = 0
        return length
    
    @staticmethod
    def load_free():
        '''
        Return the number of the next free ID
        '''
        return Momentum.load_number() + 1
    
    @staticmethod
    def load(idx):
        '''
        Load a momentum from the storage
        
        Parameters
        ----------
        idx : int
            index of the momentum in the database 'idx' > 0
        
        Returns
        -------
        momentum : Momentum
            the momentum
        '''
        
        storage = Momentum.storage
        
        #TODO: Check, for some reason some idx are given as numpy.in32 and netcdf4 is not compatible with indices given in this format!!!!!
        idx = int(idx)
        
        v = storage.ncfile.variables['momentumvelocities'][idx,:,:].astype(np.float32).copy()
        velocities = Quantity(v, nanometers / picoseconds)              
        T = storage.ncfile.variables['momentum_kinetic'][idx]
        kinetic_energy = Quantity(T, kilojoules_per_mole)
    
        momentum = Momentum(velocities=velocities, kinetic_energy=kinetic_energy)
        momentum.idx = idx

        return momentum

    @staticmethod
    def velocities_as_numpy(self, frame_indices=None, atom_indices=None):
        
        if frame_indices is None:
            frame_indices = slice(None)
        
        if atom_indices is None:
            atom_indices = slice(None)

        return self.ncfile.variables['momentumvelocities'][frame_indices,atom_indices,:].astype(np.float32).copy()
    
    @staticmethod
    def get(indices):
        return [Momentum.load(idx) for idx in range(0,Momentum.load_number())[indices] ]
            
        
    
    @staticmethod
    def _restore_netcdf(storage):
        """
        Fill in missing part after the storage has been loaded from a file and is not initialize freshly
        
        """
                
        Momentum.storage = storage
                
        return
        
    @staticmethod
    def _init_netcdf(storage):
        '''
        Initializes the associated storage to save momentums in it
        '''           
        # save associated storage in class variable for all Momentum instances to access
        
#        ncgrp = storage.ncfile.createGroup('momentum')
        
        Momentum.storage = storage
        ncgrp = storage.ncfile
        
        system = Snapshot.simulator.simulation.system
        
        # define dimensions used in momentums
        ncgrp.createDimension('momentum', 0)                       # unlimited number of momentums
        if 'atom' not in ncgrp.dimensions:
            ncgrp.createDimension('atom', system.getNumParticles())    # number of atoms in the simulated system
        
        if 'spatial' not in ncgrp.dimensions:
            ncgrp.createDimension('spatial', 3) # number of spatial dimensions        

        # define variables for momentums
        ncvar_momentumvelocities           = ncgrp.createVariable('momentumvelocities',  'f', ('momentum','atom','spatial'))
        ncvar_momentum_kinetic              = ncgrp.createVariable('momentum_kinetic',     'f', ('momentum'))

        # Define units for momentum variables.
        setattr(ncvar_momentumvelocities,  'units', 'nm/ps')
        setattr(ncvar_momentum_kinetic,     'units', 'kJ/mol')
        
        # Define long (human-readable) names for variables.
        setattr(ncvar_momentumvelocities,    "long_name", "velocities[momentum][atom][coordinate] are velocities of atom 'atom' in dimension 'coordinate' of momentum 'momentum'.")

class Snapshot(object):
    """
    Simulation snapshot.

    """
    
    # Class variables to store the global storage and the system context describing the system to be safed as snapshots
    storage = None
    simulator = None
    
    def __init__(self, context=None, coordinates=None, velocities=None, box_vectors=None, potential_energy=None, kinetic_energy=None):
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
        self.configuration = Configuration()
        self.momentum = Momentum()
        self.reversed = False

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
        return len(self.coordinates.shape[0])  

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
        this = Snapshot()
        this.momentum = self.momentum
        this.configuration = self.configuration
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
    
    def save(self, idx_configuration = None, idx_momentum = None):
        """
        Save positions, velocities, boxvectors and energies of current iteration to NetCDF file.
        
        Notes
        -----        
        We need to allow for reversed snapshots to save memory. Would be nice        
        """
        
        self.configuration.save(idx_configuration)
        self.momentum.save(idx_momentum)
        
    @staticmethod
    def load(idx_configuration = None, idx_momentum = None, reversed = False):
        '''
        Load a snapshot from the storage
        
        Parameters
        ----------
        idx : int
            index of the snapshot in the database 'idx' > 0
        
        Returns
        -------
        snapshot : Snapshot
            the snapshot
        '''
        
        #TODO: Check, for some reason some idx are given as numpy.in32 and netcdf4 is not compatible with indices given in this format!!!!!

        snapshot = Snapshot()
        snapshot.reversed = reversed
        
        if idx_configuration is not None:
            idx_c = int(idx_configuration)
            snapshot.configuration = Configuration.load(idx_c)
            
        if idx_momentum is not None:
            idx_m = int(idx_momentum)
            snapshot.momentum = Momentum.load(idx_m)

        return snapshot

    @staticmethod
    def coordinates_as_numpy(self, frame_indices=None, atom_indices=None):
        return Configuration.coordinates_as_numpy(self, frame_indices, atom_indices)
    
    @staticmethod
    def _restore_netcdf(storage):
        """
        Fill in missing part after the storage has been loaded from a file and is not initialize freshly
        
        """
                
        Snapshot.storage = storage
        Configuration.storage = storage
        Momentum.storage = storage
                
        return
        
    @staticmethod
    def _init_netcdf(storage):
        '''
        Initializes the associated storage to save snapshots in it
        '''           
        # save associated storage in class variable for all Snapshot instances to access
        # define dimensions used in snapshots

        system = Snapshot.simulator.simulation.system
        ncgrp = storage.ncfile        
        ncgrp.createDimension('atom', system.getNumParticles())    # number of atoms in the simulated system

        Configuration._init_netcdf(storage)
        Momentum._init_netcdf(storage)        
        Snapshot.storage = storage
