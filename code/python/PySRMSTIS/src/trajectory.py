'''
@author: JD Chodera
@author: JH Prinz
'''

import copy
import numpy as np

import mdtraj as md

from simtk.unit import nanosecond, picosecond, nanometers, nanometer, picoseconds, femtoseconds, femtosecond, angstroms

#=============================================================================================
# SIMULATION TRAJECTORY
#=============================================================================================

class Trajectory(list):
    """
    Simulation trajectory. Essentially a python list of snapshots

    """
    
    storage = None
    simulator = None
    
    def __init__(self, trajectory=None):
        """
        Create a simulation trajectory object

        Parameters
        ----------

        trajectory : Trajectory
            if specified, make a deep copy of specified trajectory
        """

        # Initialize list.
        list.__init__(self)
        
        self.use_lazy = True    # We assume that snapshots are immutable. That should safe a lot of time to copy trajectories

        if trajectory is not None:
            # Try to make a copy out of whatever container we were provided
            if hasattr(trajectory, 'atom_indices'):
                self.atom_indices = trajectory.atom_indices
            else:
                self.atom_indices = None
                
            if (self.use_lazy):
                self.extend(trajectory)
            else:
                for snapshot in trajectory:
                    snapshot_copy = copy.deepcopy(snapshot)                    
                    self.forward(snapshot_copy)
        else:
            self.atom_indices = None
        return

    @property
    def reversed(self):
        '''
        Returns a reversed (shallow) copy of the trajectory itself.
        '''

        t = Trajectory(self)
        t.reverse()
        return t
    
    def reverse(self):
        """
        Reverse the trajectory.

        Notes
        -----        
        We cannot handle the velocities correctly when reversing the trajectory, so velocities will no longer be meaningful.
        Kinetic energies are correctly updated, however, and path actions should be accurate.

        """
        # Reverse the order of snapshots within the trajectory.
        list.reverse(self)

        # Determine number of snapshots.
#        nsnapshots = self.__len__()
        
        # Recalculate kinetic energies for the *beginning* of each trajectory segment.
        # This makes use of the fact that the energy is (approximately) conserved over each trajectory segment, in between velocity randomizations.
        # Note that this may be a poor approximation in some cases.
#        for t in range(nsnapshots-1):
#            self[t].kinetic_energy = self[t+1].total_energy - self[t].potential_energy

        # TODO: Snapshots are immutable. Find a way to store that a snapshot has reversed velocities
        # No use reversing momenta, since we can't determine what appropriate reversed momenta should be.
        
        # We could easily indicate reversed momenta by using a minus sign in front of the index
        # this keeps everything the same and we do not need to resave the snapshots and a -idx just means take snapshot idx but invert momenta

        for frame in self:
            frame.reversed = not frame.reversed
        
        return
    
    def coordinates(self):
        """
        Return all coordinates as a numpy array
        
        Returns
        -------        
        coordinates (numpy.array(n_frames, n_atoms, 3) - numpy.array of coordinates of size number of frames 'n_frames' x number of atoms 'n_atoms' x 3 in x,y,z
        """
        
        # Make sure snapshots are stored and have an index and then add the snapshot index to the trajectory

        n_frames = self.frames     
        n_atoms = self.atoms
            
        output = np.zeros([n_frames, n_atoms, 3], np.float32)
        
        for frame_index in range(n_frames):      
            if self.atom_indices is None:
                output[frame_index,:,:] = self[frame_index].coordinates
            else:
                output[frame_index,:,:] = self[frame_index].coordinates[self.atom_indices,:]
        
        return output
    
    @property
    def frames(self):
        """
        Return the number of frames in the trajectory
        
        Returns
        -------        
        length (int) - the number of frames in the trajectory
        
        """

        return len(self)
        
    def configurations(self):
        """
        Return a list of the snapshot IDs in the trajectory
        
        Returns
        -------        
        indices (list of int) - the list of indices
        
        Notes
        -----        
        The IDs are only non-zero if the snapshots have been saved before!
        
        """
        return [f.configuration.idx for f in self]

    def momenta(self):
        """
        Return a list of the snapshot IDs in the trajectory
        
        Returns
        -------        
        indices (list of int) - the list of indices
        
        Notes
        -----        
        The IDs are only non-zero if the snapshots have been saved before!
        
        """
        return [f.momenta.idx for f in self]

    
    @property
    def atoms(self):
        """
        Return the number of atoms in the trajectory in the current view. 
        
        Returns
        -------        
        n_atoms (int) - number of atoms

        Notes
        -----        
        If a trajectory has been subsetted then this returns only the number of the view otherwise if equals the number of atoms in the snapshots stored
        
        """

        if self.atom_indices is None:
            n_atoms = self[0].coordinates.shape[0]
        else:
            n_atoms = len(self.atom_indices)
        return n_atoms    
        
    #=============================================================================================
    # LIST INHERITANCE FUNCTIONS
    #=============================================================================================
    
    def as_list(self):
        """
        Return the contained list of snapshots as a python list object
        
        Returns
        -------        
        trajectory : list of Snapshot 
            the indexed trajectory
        
        """
        return list(self)
            
    def __getslice__(self, *args, **kwargs):
        ret =  super(Trajectory, self).__getslice__(*args, **kwargs)
        if isinstance(ret, list):
            ret = Trajectory(ret)
            ret.atom_indices = self.atom_indices
            
        return ret
        
    def __getitem__(self, index):
        # Allow for numpy style of selecting several indices using a list as index parameter
        if type(index) is list:
            ret = [ super(Trajectory, self).__getitem__(i) for i in index ]
        else:
            ret = super(Trajectory, self).__getitem__(index)
                
        if isinstance(ret, list):
            ret = Trajectory(ret)
            ret.atom_indices = self.atom_indices

        return ret
    
    def __add__(self, other):        
        t = Trajectory(self)
        t.extend(other)
        return t
    
    #=============================================================================================
    # PATH ENSEMBLE FUNCTIONS
    #=============================================================================================
    
    def pathHamiltonian(self):
        """
        Compute the generalized path Hamiltonian of the trajectory.

        Parameters
        ----------
        trajectory (Trajectory) - the trajectory

        Returns
        -------        
        H : simtk.unit.Quantity with units of energy
            the generalized path Hamiltonian

        References
        ----------       
        For a description of the path Hamiltonian, see [1]:

        [1] Chodera JD, Swope WC, Noe F, Prinz JH, Shirts MR, and Pande VS. Dynamical reweighting:
        Improved estimates of dynamical properties from simulations at multiple temperatures.    
        """

        nsnapshots = len(self)
        if nsnapshots > 0:
            H = self[0].total_energy
            for snapshot_index in range(1, nsnapshots-1):
                H += self[snapshot_index].kinetic_energy
        else:
            H = 0

        return H

    def computeActivity(self):
        """
        Compute the (timeless!) activity of a given trajectory, defined in Ref. [1] as

        K[x(t)] / delta_t = delta_t \sum_{t=0}^{t_obs} \sum_{j=1}^N [r_j(t+delta_t) - r_j(t)]^2 / delta_t

        RETURNS

        K (simtk.unit) - activity K[x(t)] for the specified trajectory
        
        NOTES
        
        Can we avoid dividing and multipying by nanometers to speed up?

        """

        # Determine number of frames in trajectory.
        nframes = len(self)

        # Compute activity of component A.
        K = 0.0 * nanometers**2
        for frame_index in range(nframes-1):
            # Compute displacement of all atoms.
            delta_r = self[frame_index+1].coordinates - self[frame_index].coordinates
            # Compute contribution to activity K.
            K += ((delta_r[0:self.N,:] / nanometers)**2).sum() * (nanometers**2)

        return K 
    
    def logEquilibriumTrajectoryProbability(self):
        """
        Compute the (temperatureless!) log equilibrium probability (up to an unknown additive constant) of an unbiased trajectory evolved according to Verlet dynamics with Andersen thermostatting.

        ARGUMENTS
        trajectory (Trajectory) - the trajectory

        Returns
        -------        
        log_q : float
            the log equilibrium probability of the trajectory divided by the inverse temperature beta
        
        NOTES
        This might be better places into the trajectory class. The trajectory should know the system and ensemble? and so it is not necessarily 
        TPS specific

        """

        nsnapshots = len(self)
        log_q = - self[0].total_energy
        for snapshot_index in range(1, nsnapshots-1):
            log_q += - self[snapshot_index].kinetic_energy

        return log_q

    
    #=============================================================================================
    # UTILITY FUNCTIONS
    #=============================================================================================

    def subset(self, atom_indices):
        """
        Reduce the view of the trajectory to a subset of atoms specified. This is only a view, no data will be changed or copied.
        
        Returns
        -------        
        trajectory : Trajectory
            the trajectory showing the subsets of atoms
        """        

        t = Trajectory(self)
        t.atom_indices = atom_indices
        return t

    @property
    def solute(self):
        """
        Reduce the view of the trajectory to a subset of solute atoms specified in the associated Simulator
        
        Returns
        -------        
        trajectory : Trajectory
            the trajectory showing the subsets of solute atoms
        """        

        return self.subset(Trajectory.simulator.solute_indices)

    def md(self):
        """
        Construct a mdtraj.Trajectory object from the Trajectory itself
        
        Returns
        -------        
        trajectory : mdtraj.Trajectory
            the trajectory
        """        

        output = self.coordinates()
        topology = self.md_topology()
                                                 
        return md.Trajectory(output, topology)             
                
    
    def md_topology(self):
        """
        Return a mdtraj.Topology object representing the topology of the current view of the trajectory
        
        Returns
        -------        
        topology : mdtraj.Topology
            the topology
        """        

        topology = md.Topology.from_openmm(Trajectory.simulator.simulation.topology)
        
        if self.atom_indices is not None:
            topology = topology.subset(self.atom_indices)       
        
        return topology

    #=============================================================================================
    # STORAGE FUNCTIONS
    #=============================================================================================
    
    def save(self):
        """
        Add the current state of the trajectory in the database. If nothing has changed then the trajectory gets stored using the same snapshots as before. Saving lots of diskspace
        
        """
        
        ncfile = Trajectory.storage.ncfile
        idx = Trajectory.load_free( )

        # Make sure snapshots are stored and have an index and then add the snapshot index to the trajectory

        nframes = len(self)
        for frame_index in range(nframes):      
            frame = self[frame_index]         
            frame.save()
            ncfile.variables['trajectory_configuration_idx'][idx,frame_index] = frame.configuration.idx         
            ncfile.variables['trajectory_momentum_idx'][idx,frame_index] = frame.momentum.idx                     
            ncfile.variables['trajectory_momentum_reversed'][idx,frame_index] = frame.reversed                     
             
        ncfile.variables['trajectory_length'][idx] = nframes
                
        return 

    @staticmethod
    def load_momentum_indices(idx):
        '''
        Load trajectory indices for trajectory with ID 'idx' from the storage
        
        ARGUMENTS
        
        idx (int) - ID of the trajectory
        
        Returns
        -------        
        trajectory (list of int) - trajectory indices
        '''    
        
        
        
        length = Trajectory.load_length(idx)
        return Trajectory.storage.ncfile.variables['trajectory_momentum_idx'][idx,:length]    

    
    @staticmethod
    def load_configuration_indices(idx):
        '''
        Load trajectory indices for trajectory with ID 'idx' from the storage
        
        ARGUMENTS
        
        idx (int) - ID of the trajectory
        
        Returns
        -------        
        trajectory (list of int) - trajectory indices
        '''    
        
        
        
        length = Trajectory.load_length(idx)
        return Trajectory.storage.ncfile.variables['trajectory_configuration_idx'][idx,:length]    

    @staticmethod
    def load_indices(idx):
        '''
        Load trajectory indices for trajectory with ID 'idx' from the storage
        
        ARGUMENTS
        
        idx (int) - ID of the trajectory
        
        Returns
        -------        
        trajectory (list of int) - trajectory indices
        '''    
        
        length = Trajectory.load_length(idx)
        return (Trajectory.storage.ncfile.variables['trajectory_configuration_idx'][idx,:length],
                   Trajectory.storage.ncfile.variables['trajectory_momentum_idx'][idx,:length])

    @staticmethod
    def load_length(idx):
        '''
        Return the length of a trajectory from the storage
        
        Parameters
        ----------
        idx : int
            index of the trajectory
            
        Returns
        -------
        length : int
            Number of frames in the trajectory
        '''
        return Trajectory.storage.ncfile.variables['trajectory_length'][idx]
    
    @staticmethod
    def load(idx):
        '''
        Return a trajectory from the storage
        
        Parameters
        ----------
        idx : int
            index of the trajectory
            
        Returns
        -------
        trajectory : Trajectory
            the trajectory
        '''        
        frames_c, frames_m = Trajectory.load_indices(idx)
        return Trajectory.from_indices(frames_c, frames_m)

    @staticmethod
    def from_indices(frames_configuration, frames_momentum):
        '''
        Return a trajectory from the storage constructed from a list of snapshot indices
        
        Parameters
        ----------
        frames : list of int
            list of snapshot indices to be used to generate the trajectory
            
        Returns
        -------
        trajectory : Trajectory
            the trajectory
        '''
        trajectory = Trajectory()

        for frame in zip(frames_configuration, frames_momentum):                
            snapshot = Trajectory.storage.snapshot( frame[0], frame[1] )
            trajectory.append(snapshot)
        
        return trajectory
    
    @staticmethod
    def load_number():
        '''
        Return the number of trajectories in the storage
        
        Returns
        -------
        number : int
            number of trajectories in the storage. Their indexing starts with 1.
        '''
        length = int(len(Trajectory.storage.ncfile.dimensions['trajectory'])) - 1
        if length < 0:
            length = 0
        return length
    
    @staticmethod
    def load_free():
        '''
        Return the number of the next free ID
        
        Returns
        -------
        index : int
            the number of the next free index in the storage. Used to store a new snapshot.
        '''
        return Trajectory.load_number() + 1

    @staticmethod
    def load_all_momentum_indices():
        '''
        Return a list of frame indices for all trajectories in the storage
        
        Returns
        -------        
        list : list of list of int
            a list of list of frame IDs for all trajectories.
        '''
        ind = Trajectory.storage.ncfile.variables['trajectory_momentum_idx'][:,:].astype(np.int32).copy()
        lengths = Trajectory.storage.ncfile.variables['trajectory_length'][:].astype(np.int32).copy()
        n_traj = Trajectory.load_number()
        
        return [ ind[i,:lengths[i]] for i in range(1,n_traj + 1) ]    
    
    @staticmethod
    def load_all_configuration_indices():
        '''
        Return a list of frame indices for all trajectories in the storage
        
        Returns
        -------        
        list : list of list of int
            a list of list of frame IDs for all trajectories.
        '''
        ind = Trajectory.storage.ncfile.variables['trajectory_configuration_idx'][:,:].astype(np.int32).copy()
        lengths = Trajectory.storage.ncfile.variables['trajectory_length'][:].astype(np.int32).copy()
        n_traj = Trajectory.load_number()
        
        return [ ind[i,:lengths[i]] for i in range(1,n_traj + 1) ]            
    
    @staticmethod
    def _restore_netcdf(storage):
        """
        Fill in missing part after the storage has been loaded from a file and is not initialize freshly
        
        """
        Trajectory.storage = storage        
    
    @staticmethod
    def _init_netcdf(storage):        
        """
        Initialize the associated storage to allow for trajectory storage
        
        """        

        # save associated storage in class variable for all Trajectory instances to access
        Trajectory.storage = storage
        ncfile = storage.ncfile
                
        nframes = Trajectory.simulator.n_frames_max
        
        # define dimensions used in trajectories
        ncfile.createDimension('trajectory', 0)                 # unlimited number of iterations
        ncfile.createDimension('max_frames', nframes)     # number of maximal frames per trajectory

        # Create variables for trajectories        
        ncvar_trajectory_configuration_idx  = ncfile.createVariable('trajectory_configuration_idx', 'u4', ('trajectory','max_frames'))
        ncvar_trajectory_momentum_idx       = ncfile.createVariable('trajectory_momentum_idx', 'u4', ('trajectory','max_frames'))
        ncvar_trajectory_momentum_reversed  = ncfile.createVariable('trajectory_momentum_reversed', 'b', ('trajectory', 'max_frames'))
        ncvar_trajectory_path_hamiltonian   = ncfile.createVariable('trajectory_path_hamiltonians', 'f', ('trajectory'))
        ncvar_trajectory_length             = ncfile.createVariable('trajectory_length', 'u4', ('trajectory'))

        # Define units for snapshot variables.
        setattr(ncvar_trajectory_path_hamiltonian,      'units', 'none')
        setattr(ncvar_trajectory_configuration_idx,     'units', 'none')
        setattr(ncvar_trajectory_momentum_idx,          'units', 'none')
        setattr(ncvar_trajectory_length,                'units', 'none')
        setattr(ncvar_trajectory_momentum_reversed,     'units', 'none')
        
        # Define long (human-readable) names for variables.
        setattr(ncvar_trajectory_configuration_idx,    "long_name", "trajectory[trajectory][frame] is the snapshot index (0..nspanshots-1) of frame 'frame' of trajectory 'trajectory'.")
        setattr(ncvar_trajectory_momentum_idx,         "long_name", "trajectory[trajectory][frame] is the snapshot index (0..nspanshots-1) of frame 'frame' of trajectory 'trajectory'.")