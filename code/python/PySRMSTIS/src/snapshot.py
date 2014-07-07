'''

@author: JD Chodera
'''

import copy
import numpy

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
    context = None
    
    def __init__(self, context=None, coordinates=None, velocities=None, box_vectors=None, potential_energy=None, kinetic_energy=None):
        """
        Create a simulation snapshot from either an OpenMM context or individually-specified components.

        OPTIONAL ARGUMENTS

        context (simtk.chem.openContext) - if not None, the current state will be queried to populate simulation snapshot; otherwise, can specify individual components (default: None)
        coordinates (simtk.unit.Quantity wrapping Nx3 numpy array of dimension length) - atomic coordinates (default: None)
        velocities (simtk.unit.Quantity wrapping Nx3 numpy array of dimension length) - atomic velocities (default: None)
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
        if numpy.any(numpy.isnan(self.coordinates)):
            raise Exception("Some coordinates became 'nan'; simulation is unstable or buggy.")

        return

    @property
    def total_energy(self):
        return self.kinetic_energy + self.potential_energy
    
    def save(self):
        """
        Save positions, states, and energies of current iteration to NetCDF file.
        
        """
        
        if self.idx > 0:
            idx = Snapshot.storage.snapshot_idx
    
            # Store snapshot.
            self.ncfile.variables['snapshot_coordinates'][idx,:,:] = (self.coordinates / nanometers).astype(numpy.float32)
            self.ncfile.variables['snapshot_velocities'][idx,:,:] = (self.velocities / (nanometers / picoseconds)).astype(numpy.float32)
            self.ncfile.variables['snapshot_potential'][idx] = self.potential_energy / kilojoules_per_mole                                
            self.ncfile.variables['snapshot_kinetic'][idx] = self.kinetic_energy / kilojoules_per_mole
            
            # increase snapshout counter and store for later reference in Snapshot object
            
            self.idx = idx
            Snapshot.storage.snapshot_idx += 1        
    
            # Force sync to disk to avoid data loss.
            self.ncfile.sync()

        return
    
    @staticmethod
    def _init_netcdf(storage):
        
        # save associated storage in class variable for all Snapshot instances to access
        Snapshot.storage = storage
        ncfile = storage.ncfile
        
        storage.snapshot_idx = 1;
        
        system = Snapshot.context.getSystem()
        
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