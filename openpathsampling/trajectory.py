'''
@author: JD Chodera
@author: JH Prinz
'''

import copy

import numpy as np
import mdtraj as md
import simtk.unit as u

import openpathsampling as paths

#=============================================================================================
# SIMULATION TRAJECTORY
#=============================================================================================



class Trajectory(list):
    """
    Simulation trajectory. Essentially a python list of snapshots

    """
    
    engine = None
    use_lazy = True    # We assume that snapshots are immutable. That should safe a lot of time to copy trajectories


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

        self.path_probability = None # For future uses

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
                    self.append(snapshot_copy)
        else:
            self.atom_indices = None

    def __str__(self):
        return 'Trajectory[' + str(len(self)) + ']'

    def __repr__(self):
        return 'Trajectory[' + str(len(self)) + ']'

    @property
    def reversed(self):
        '''
        Returns a reversed (shallow) copy of the trajectory itself. Effectively
        creates a new Trajectory object and then fills it with shallow reversed copies
        of the contained snapshots.

        Returns
        -------
        Trajectory()
            the reversed trajectory
        '''

        return Trajectory([snap.reversed_copy() for snap in reversed(self)])

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
        n_spatial = self.spatial
            
        output = np.zeros([n_frames, n_atoms, n_spatial], np.float32)
        
        for frame_index in range(n_frames):      
            if self.atom_indices is None:
                output[frame_index,:,:] = self[frame_index].coordinates
            else:
                output[frame_index,:,:] = self[frame_index].coordinates[self.atom_indices,:]

        return output
    
    @property
    def frames(self):
        """
        Return the number of frames in the trajectory.
        
        Returns
        -------        
        length (int) - the number of frames in the trajectory

        Notes
        -----
        Might be removed in later versions len(trajectory) is more intuitive
        
        """

        return len(self)

    def configurations(self):
        """
        Return a list of the snapshots in the trajectory

        Returns
        -------
        list of Configuration
            the list of Configuration objects

        """
        return [f.configuration for f in self]


    def momenta(self):
        """
        Return a list of the Momentum objects in the trajectory

        Returns
        -------
        list of Momentum()
            the list of Momentum objects
        """
        return [f.momenta for f in self]

    @property
    def spatial(self):
        if self.topology is None:
            n_spatial = self[0].coordinates.shape[1]
        else:
            n_spatial = self.topology.n_spatial

        return n_spatial

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

    def __getslice__(self, *args, **kwargs):
        ret =  list.__getslice__(self, *args, **kwargs)
        if isinstance(ret, list):
            ret = Trajectory(ret)
            ret.atom_indices = self.atom_indices
            
        return ret
        
    def __getitem__(self, index):
        # Allow for numpy style of selecting several indices using a list as index parameter
        if type(index) is list:
            ret = [ list.__getitem__(self, i) for i in index ]
        else:
            ret = list.__getitem__(self, index)
                
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
        K = 0.0 * u.nanometers**2
        for frame_index in range(nframes-1):
            # Compute displacement of all atoms.
            delta_r = self[frame_index+1].coordinates - self[frame_index].coordinates
            # Compute contribution to activity K.
            K += ((delta_r[0:self.N,:] / u.nanometers)**2).sum() * (u.nanometers**2)

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
        Reduce the view of the trajectory to a subset of solute atoms
        specified in the associated DynamicsEngine
        
        Returns
        -------        
        trajectory : Trajectory
            the trajectory showing the subsets of solute atoms
        """

        #TODO: To remove the dependency of the dynamics engine we need to get the information
        #TODO: about the solute_indices from somewhere else, preferrably the topology?

        if Trajectory.engine is None:
            raise ValueError("No engine specified to get solute_indices from !")

        return self.subset(Trajectory.engine.solute_indices)

    def full(self):
        """
        Return a view of the trajectory with all atoms

        Returns
        -------
        trajectory : Trajectory
            the trajectory showing the subsets of solute atoms
        """
        return self.subset(None)

    def md(self, topology = None):
        """
        Construct a mdtraj.Trajectory object from the Trajectory itself

        Parameters
        ----------
        topology : mdtraj.Topology()
            If not None this topology will be used to construct the mdtraj
            objects otherwise the topology object will be taken from the
            configurations in the trajectory snapshots.
        
        Returns
        -------        
        trajectory : mdtraj.Trajectory
            the trajectory
        """

        if topology is None:
            topology = self.md_topology()

        output = self.coordinates()

        return md.Trajectory(output, topology)

    @staticmethod
    def from_mdtraj(mdtrajectory):
        """
        Construct a Trajectory object from an mdtraj.Trajectory object

        Parameters
        ----------
        mdtrajectory : mdtraj.Trajectory
            Input mdtraj.Trajectory

        Returns
        -------
        Trajectory
        """
        trajectory = Trajectory()
        empty_momentum = paths.Momentum()
        empty_momentum.velocities = None
        for frame_num in range(mdtrajectory.n_frames):
            # mdtraj trajectories only have coordinates and box_vectors
            coord = u.Quantity(mdtrajectory.xyz[frame_num], u.nanometers)
            if mdtrajectory.unitcell_vectors is not None:
                box_v = u.Quantity(mdtrajectory.unitcell_vectors[frame_num],
                                 u.nanometers)
            else:
                box_v = None
            config = paths.Configuration(coordinates=coord, box_vectors=box_v)

            snap = paths.Snapshot(configuration=config, momentum=empty_momentum)
            trajectory.append(snap)

        return trajectory

    @property
    def topology(self):
        """
        Return a Topology object representing the topology of the
        current view of the trajectory

        Returns
        -------
        topology : opentis.Topology
            the topology object

        Notes
        -----
        This is taken from the configuration of the first frame.
        """

        topology = None

        if len(self) > 0 and self[0].topology is not None:
            # if no topology is defined
            topology = self[0].topology

            if self.atom_indices is not None:
                topology = topology.subset(self.atom_indices)

        return topology

    def md_topology(self):
        """
        Return a mdtraj.Topology object representing the topology of the
        current view of the trajectory
        
        Returns
        -------        
        topology : mdtraj.Topology
            the topology

        Notes
        -----
        This is taken from the configuration of the first frame.
        Use topology.md instead
        TODO: Should be removed
        """        

        return self.topology.md
