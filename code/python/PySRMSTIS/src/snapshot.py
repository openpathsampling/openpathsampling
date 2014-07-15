'''

@author: JD Chodera
'''

import copy
import numpy as np
import mdtraj as md

from simtk.unit import nanosecond, picosecond, nanometers, nanometer, picoseconds, femtoseconds, femtosecond, kilojoules_per_mole, Quantity


#=============================================================================================
# SIMULATION SNAPSHOT 
#=============================================================================================

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

        OPTIONAL ARGUMENTS

        context (simtk.chem.openContext) - if not None, the current state will be queried to populate simulation snapshot; otherwise, can specify individual components (default: None)
        coordinates (simtk.unit.Quantity wrapping Nx3 np array of dimension length) - atomic coordinates (default: None)
        velocities (simtk.unit.Quantity wrapping Nx3 np array of dimension length) - atomic velocities (default: None)
        box_vectors - periodic box vectors (default: None)
        potential_energy (simtk.unit.Quantity of units energy/mole) - potential energy at current timestep (default: None)
        kinetic_energy (simtk.unit.Quantity of units energy/mole) - kinetic energy at current timestep (default: None)
        
        """
        
        self.idx = 0        # potential idx in a netcdf storage, if 0 then not stored yet. Attention! Cannot be stored in 2 repositories at the same time

        if context is not None:
            # Get current state from OpenMM Context object.
            state = context.getState(getPositions=True, getVelocities=True, getEnergy=True)
            
            # Store the associated context
            self.context = context
            
            # Populate current snapshot data.
            self.coordinates = state.getPositions(asNumpy=True)
            self.velocities = state.getVelocities(asNumpy=True)
            self.box_vectors = state.getPeriodicBoxVectors() # TODO: set asNumpy=True once bug in OpenMM is fixed
            self.potential_energy = state.getPotentialEnergy()
            self.kinetic_energy = state.getKineticEnergy()
        else:
            if coordinates is not None: self.coordinates = copy.deepcopy(coordinates)
            if velocities is not None: self.velocities = copy.deepcopy(velocities)
            if box_vectors is not None: self.box_vectors = copy.deepcopy(box_vectors)
            if potential_energy is not None: self.potential_energy = copy.deepcopy(potential_energy)
            if kinetic_energy is not None: self.kinetic_energy = copy.deepcopy(kinetic_energy)                       

        # Check for nans in coordinates, and raise an exception if something is wrong.
        if np.any(np.isnan(self.coordinates)):
            raise Exception("Some coordinates became 'nan'; simulation is unstable or buggy.")

        return

    @property
    def atoms(self):
        '''
        returns the number of atoms in the snapshot
        '''   
        return len(self.coordinates.shape[0])  

    @property
    def total_energy(self):
        '''
        returns the total energy (sum of potential and kinetic) of the snapshot
        '''   
        return self.kinetic_energy + self.potential_energy
    
    #=============================================================================================
    # Utility functions
    #=============================================================================================

    def md(self):
        '''
        returns a mdtraj Trajectory object that contains only one frame
        
        NOTES
        
        Rather slow since the topology has to be made each time. Try to avoid it
        '''        
        
        n_atoms = self.atoms
                            
        output = np.zeros([1, n_atoms, 3], np.float32)
        output[0,:,:] = self.coordinates
        
        topology = self.md_topology()
                                                         
        return md.Trajectory(output, topology)      

    
    def md_topology(self): 
        '''
        returns a mdtraj topology object that can be used with the stored snapshot
        '''   
        return md.Topology.from_openmm(Snapshot.simulator.simulation.topology)
 
    #=============================================================================================
    # Storage functions
    #=============================================================================================
    
    def save(self, idx = None):
        """
        Save positions, states, and energies of current iteration to NetCDF file.
        
        NOTES
        
        We need to allow for reversed snapshots to save memory. Would be nice
        
        """
        
        storage = Snapshot.storage
        
        if self.idx == 0:
            if idx is None:
                idx = storage.snapshot_idx
                storage.snapshot_idx += 1        
    
            # Store snapshot.
            storage.ncfile.variables['snapshot_coordinates'][idx,:,:] = (self.coordinates / nanometers).astype(np.float32)
            storage.ncfile.variables['snapshot_velocities'][idx,:,:] = (self.velocities / (nanometers / picoseconds)).astype(np.float32)
            storage.ncfile.variables['snapshot_potential'][idx] = self.potential_energy / kilojoules_per_mole                                
            storage.ncfile.variables['snapshot_kinetic'][idx] = self.kinetic_energy / kilojoules_per_mole
            
            # increase snapshout counter and store for later reference in Snapshot object
            
            self.idx = idx
    
            # Force sync to disk to avoid data loss.
            storage.ncfile.sync()

        return

    @staticmethod
    def load_number():
        '''
        Load a snapshot from the storage
        
        RETURNS
        
        number (int) - number of stored snapshots
        '''
        length = len(Snapshot.storage.ncfile.dimensions['snapshot']) - 1
        if length < 0:
            length = 0
        return length
        
    
    @staticmethod
    def load(idx):
        '''
        Load a snapshot from the storage
        
        ARGUMENTS
        
        idx (int) - index of the snapshot in the database 'idx' > 0
        
        RETURNS
        
        snapshot (Snapshot) - the snapshot
        '''
        
        storage = Snapshot.storage
        
#        print type(idx)
        #TODO: Check, for some reason some idx are given as numpy.in32 and netcdf4 is not compatible with indices given in this format!!!!!
        idx = int(idx)
        
        x = storage.ncfile.variables['snapshot_coordinates'][idx,:,:].astype(np.float32).copy()
        coordinates = Quantity(x, nanometers)                
        v = storage.ncfile.variables['snapshot_velocities'][idx,:,:].astype(np.float32).copy()
        velocities = Quantity(v, nanometers / picoseconds)                
        V = storage.ncfile.variables['snapshot_potential'][idx]
        potential_energy = Quantity(V, kilojoules_per_mole)
        T = storage.ncfile.variables['snapshot_kinetic'][idx]
        kinetic_energy = Quantity(T, kilojoules_per_mole)
        snapshot = Snapshot(coordinates=coordinates, velocities=velocities, kinetic_energy=kinetic_energy, potential_energy=potential_energy)
        snapshot.idx = idx
        return snapshot

    @staticmethod
    def coordinates_as_numpy(self, frame_indices=None, atom_indices=None):
        
        if frame_indices is None:
            frame_indices = slice(None)
        
        if atom_indices is None:
            atom_indices = slice(None)

        return self.ncfile.variables['snapshot_coordinates'][frame_indices,atom_indices,:].astype(np.float32).copy()
    
    @staticmethod
    def _restore_netcdf(storage):
        """
        Fill in missing part after the storage has been loaded from a file and is not initialize freshly
        
        NOTES
        
        We need to allow for reversed snapshots to save memory. Would be nice
        
        """
        
        Snapshot.storage = storage
        
        storage.snapshot_idx = int(storage.ncfile.variables['snapshot_coordinates'].shape[0])
        
        return
        
    @staticmethod
    def _init_netcdf(storage):
        '''
        initializes the associates storage to save snapshots in it
        '''           
        # save associated storage in class variable for all Snapshot instances to access
        Snapshot.storage = storage
        ncfile = storage.ncfile
        
        storage.snapshot_idx = 1;
        
        system = Snapshot.simulator.simulation.system
        
        # define dimensions used in snapshots
        ncfile.createDimension('snapshot', 0)                       # unlimited number of snapshots
        ncfile.createDimension('atom', system.getNumParticles())    # number of atoms in the simulated system
        ncfile.createDimension('spatial', 3)                        # number of spatial dimensions

        # define variables for snapshots
        ncvar_snapshot_coordinates          = ncfile.createVariable('snapshot_coordinates', 'f', ('snapshot','atom','spatial'))
        ncvar_snapshot_velocities           = ncfile.createVariable('snapshot_velocities',  'f', ('snapshot','atom','spatial'))
        ncvar_snapshot_potential            = ncfile.createVariable('snapshot_potential',   'f', ('snapshot'))
        ncvar_snapshot_kinetic              = ncfile.createVariable('snapshot_kinetic',     'f', ('snapshot'))

        # Define units for snapshot variables.
        setattr(ncvar_snapshot_coordinates, 'units', 'nm')
        setattr(ncvar_snapshot_velocities,  'units', 'nm/ps')
        setattr(ncvar_snapshot_potential,   'units', 'kJ/mol')
        setattr(ncvar_snapshot_kinetic,     'units', 'kJ/mol')
        
        # Define long (human-readable) names for variables.
        setattr(ncvar_snapshot_coordinates,    "long_name", "coordinates[snapshot][atom][coordinate] are coordinate of atom 'atom' in dimension 'coordinate' of snapshot 'snapshot'.")
        setattr(ncvar_snapshot_velocities,    "long_name", "velocities[snapshot][atom][coordinate] are velocities of atom 'atom' in dimension 'coordinate' of snapshot 'snapshot'.")