'''
@author: JD Chodera
@author: JH Prinz
'''

import copy
from simtk.unit import nanosecond, picosecond, nanometers, nanometer, picoseconds, femtoseconds, femtosecond, angstroms

#=============================================================================================
# SIMULATION TRAJECTORY
#=============================================================================================

class Trajectory(list):
    """
    Simulation trajectory. Essentially a python list of snapshots

    """
    
    storage = None
    nframes = 1000
    
    def __init__(self, trajectory=None):
        """
        Create a simulation trajectory object

        OPTIONAL ARGUMENTS

        trajectory (Trajectory) - if specfiied, make a deep copy of specified trajectory
        
        """

        # Initialize list.
        list.__init__(self)
        
        self.use_lazy = True    # We assume that snapshots are immutable. That should safe a lot of time to copy trajectories

        if trajectory is not None:
            # Try to make a copy out of whatever container we were provided
            for snapshot in trajectory:
                if (self.use_lazy):
                    self.append(snapshot)
                else:
                    snapshot_copy = copy.deepcopy(snapshot)                    
                    self.append(snapshot_copy)

        return

    def reverse(self):
        """
        Reverse the trajectory.

        NOTE

        We cannot handle the velocities correctly when reversing the trajectory, so velocities will no longer be meaningful.
        Kinetic energies are correctly updated, however, and path actions should be accurate.

        """
        # Reverse the order of snapshots within the trajectory.
        list.reverse(self)

        # Determine number of snapshots.
        nsnapshots = self.__len__()
        
        # Recalculate kinetic energies for the *beginning* of each trajectory segment.
        # This makes use of the fact that the energy is (approximately) conserved over each trajectory segment, in between velocity randomizations.
        # Note that this may be a poor approximation in some cases.
        for t in range(nsnapshots-1):
            self[t].kinetic_energy = self[t+1].total_energy - self[t].potential_energy

        # No use reversing momenta, since we can't determine what appropriate reversed momenta should be.
        
        return
    
    #=============================================================================================
    # PATH ENSEMBLE FUNCTIONS
    #=============================================================================================
    
    def pathHamiltonian(self):
        """
        Compute the generalized path Hamiltonian of the trajectory.

        ARGUMENTS

        trajectory (Trajectory) - the trajectory

        RETURNS

        H (simtk.unit.Quantity with units of energy) - the generalized path Hamiltonian

        REFERENCES

        For a description of the path Hamiltonian, see [1]:

        [1] Chodera JD, Swope WC, Noe F, Prinz JH, Shirts MR, and Pande VS. Dynamical reweighting:
        Improved estimates of dynamical properties from simulations at multiple temperatures.    

        """

        nsnapshots = len(self)
        H = self[0].total_energy
        for snapshot_index in range(1, nsnapshots-1):
            H += self[snapshot_index].kinetic_energy

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

        RETURNS

        log_q (float) - the log equilibrium probability of the trajectory divided by the inverse temperature beta
        
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
    
    def write_pdb(self, filename, atoms):
        """
        Write out a trajectory as a multi-model PDB files.
    
        ARGUMENTS
           filename (string) - name of PDB file to be written
           trajectory (Trajectory)
           atoms (list of dict) - parsed PDB file ATOM entries from read_pdb() or construct_atom_list(); coordinates will be modified
        """
    
        # Create file.
        outfile = open(filename, 'w')
    
        nframes = len(self)
    
        # Write trajectory as models
        for frame_index in range(nframes):
            outfile.write("MODEL     %4d\n" % (frame_index+1))
            coordinates = self[frame_index].coordinates
            # Write ATOM records.
            for (index, atom) in enumerate(atoms):
                atom["x"] = "%8.3f" % (coordinates[index,0] / angstroms)
                atom["y"] = "%8.3f" % (coordinates[index,1] / angstroms)
                atom["z"] = "%8.3f" % (coordinates[index,2] / angstroms)
                outfile.write('ATOM  %(serial)5s %(atom)4s%(altLoc)c%(resName)3s %(chainID)c%(Seqno)5s   %(x)8s%(y)8s%(z)8s\n' % atom)
    
            outfile.write("ENDMDL\n")
            
        # Close file.
        outfile.close()
    
        return
    
    def save(self):
        """
        Add the current state of the trajectory in the database. If nothing has changed then the trajectory gets stored using the same snapshots as before. Saving lots of diskspace
        
        """
        
        ncfile = Trajectory.storage.ncfile
        idx = Trajectory.storage.trajectory_idx

        # Make sure snapshots are stored and have an index and then add the snapshot index to the trajectory

        nframes = len(self)
        for frame_index in range(nframes):      
            frame = self[frame_index]         
            frame.save()
            ncfile.variables['trajectory_idx'][idx,frame_index] = frame.idx         
             
        ncfile.variables['trajectory_length'][idx] = nframes
        
        Trajectory.storage.trajectory_idx += 1
        
        return 
       
    @staticmethod
    def _init_netcdf(storage):        

        # save associated storage in class variable for all Trajectory instances to access
        Trajectory.storage = storage
        ncfile = storage.ncfile
        
        storage.trajectory_idx = 1;
        
        # define dimensions used in trajectories
        ncfile.createDimension('trajectory', 0)                 # unlimited number of iterations
        ncfile.createDimension('frame', Trajectory.nframes)     # number of maximal frames per trajectory

        # Create variables for trajectories        
        ncvar_trajectory_idx                = ncfile.createVariable('trajectory_idx', 'u4', ('trajectory','frame'))
        ncvar_trajectory_length             = ncfile.createVariable('trajectory_length', 'u4', ('trajectory'))
        ncvar_trajectory_path_hamiltonian   = ncfile.createVariable('path_hamiltonians', 'f', ('trajectory'))

        # Define units for snapshot variables.
        setattr(ncvar_trajectory_path_hamiltonian,      'units', 'none')
        setattr(ncvar_trajectory_idx,                   'units', 'none')
        setattr(ncvar_trajectory_length,                'units', 'none')
        
        # Define long (human-readable) names for variables.
        setattr(ncvar_trajectory_idx,    "long_name", "trajectory[trajectory][frame] is the snapshot index (0..nspanshots-1) of frame 'frame' of trajectory 'trajectory'.")