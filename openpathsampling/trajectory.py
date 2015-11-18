"""
@author: JD Chodera
@author: JH Prinz
"""

import numpy as np
import mdtraj as md
import simtk.unit as u

import openpathsampling as paths
from openpathsampling.netcdfplus import StorableObject


# =============================================================================================
# SIMULATION TRAJECTORY
# =============================================================================================


class Trajectory(list, StorableObject):
    """
    Simulation trajectory. Essentially a python list of snapshots
    """

    engine = None

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
        StorableObject.__init__(self)

        self.path_probability = None  # For future uses

        if trajectory is not None:
            # Try to make a copy out of whatever container we were provided
            if hasattr(trajectory, 'atom_indices'):
                self.atom_indices = trajectory.atom_indices
            else:
                self.atom_indices = None
            if type(trajectory) is Trajectory:
                self.extend(trajectory.iter_proxies())
            else:
                self.extend(trajectory)
        else:
            self.atom_indices = None

    def extend(self, iterable):
        if type(iterable) is Trajectory:
            list.extend(self, iterable.iter_proxies())
        else:
            list.extend(self, iterable)

    def __str__(self):
        return 'Trajectory[' + str(len(self)) + ']'

    def __repr__(self):
        return 'Trajectory[' + str(len(self)) + ']'

    def map(self, fnc, allow_fast=True):
        """
        This runs a function and tries to be fast.

        Fast here means that functions that are purely based on CVs can be
        evaluated without actually loading the real Snapshot object. This
        functions tries to do that and if it fails it does it the usual way
        and creates the snapshot object. This bears the possibility that
        the function uses the fake snapshots and returns a non-sense value.
        It is up to the user to make sure this will not happen.
        """

        if allow_fast:
            try:
                return [fnc(frame) for frame in list.__iter__(self)]
            except:
                pass

        return [fnc(frame) for frame in self]

    @property
    def reversed(self):
        """
        Returns a reversed (shallow) copy of the trajectory itself. Effectively
        creates a new Trajectory object and then fills it with shallow reversed copies
        of the contained snapshots.

        Returns
        -------
        Trajectory()
            the reversed trajectory
        """

        return Trajectory([snap for snap in reversed(self)])

    def prepend(self, snapshot):
        """
        Prepend a snapshot

        Just convenience method to replace insert(0, snapshot)
        """
        self.insert(0, snapshot)
        # And a generation of scientist-programmers who grew up learning
        # "OPS trajectories are just Python lists" scream in pain when they
        # find this after googling "python list.prepend not working".
        # (Blame JHP. This was his doing.)

    def coordinates(self):
        """
        Return all coordinates as a numpy array
        
        Returns
        -------        
        coordinates : numpy.ndarray((n_frames, n_atoms, 3))
            numpy.array of coordinates of size number of
            frames 'n_frames' x number of atoms 'n_atoms' x 3 in x,y,z
        """

        # Make sure snapshots are stored and have an index and then add the snapshot index to the trajectory

        n_frames = len(self)
        n_atoms = self.n_atoms
        n_spatial = self.spatial

        output = np.zeros([n_frames, n_atoms, n_spatial], np.float32)

        for frame_index in range(n_frames):
            if self.atom_indices is None:
                output[frame_index, :, :] = self[frame_index].coordinates
            else:
                output[frame_index, :, :] = self[frame_index].coordinates[self.atom_indices, :]

        return output

    def xyz(self):
        n_frames = len(self)
        n_atoms = self.n_atoms
        n_spatial = self.spatial

        output = np.zeros([n_frames, n_atoms, n_spatial], np.float32)

        for frame_index in range(n_frames):
            if self.atom_indices is None:
                output[frame_index, :, :] = self[frame_index].xyz
            else:
                output[frame_index, :, :] = self[frame_index].xyz[self.atom_indices, :]

        return output

    @property
    def n_snapshots(self):
        """
        Return the number of frames in the trajectory.
        
        Returns
        -------        
        length (int) - the number of frames in the trajectory

        Notes
        -----
        Might be removed in later versions for len(trajectory) is more pythonic

        See also
        --------
        n_frames, len

        """

        return len(self)

    @property
    def n_frames(self):
        """
        Return the number of frames in the trajectory.

        Returns
        -------
        length (int) - the number of frames in the trajectory

        See also
        --------
        n_snapshots, len

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
    def n_atoms(self):
        """
        Return the number of atoms in the trajectory in the current view. 
        
        Returns
        -------        
        n_atoms : int
            number of atoms

        Notes
        -----        
        If a trajectory has been subsetted then this returns only the number
        of the view otherwise if equals the number of atoms in the snapshots stored
        
        """

        if self.atom_indices is None:
            n_atoms = self[0].xyz.shape[0]
        else:
            n_atoms = len(self.atom_indices)
        return n_atoms

    # =============================================================================================
    # LIST INHERITANCE FUNCTIONS
    # =============================================================================================

    def __getslice__(self, *args, **kwargs):
        ret = list.__getslice__(self, *args, **kwargs)
        if type(ret) is list:
            ret = Trajectory(ret)
            ret.atom_indices = self.atom_indices

        return ret

    def __hash__(self):
        return object.__hash__(self)

    def __getitem__(self, index):
        # Allow for numpy style of selecting several indices using a list as index parameter
        if hasattr(index, '__iter__'):
            ret = [list.__getitem__(self, i) for i in index]
        else:
            ret = list.__getitem__(self, index)

        if type(ret) is list:
            ret = Trajectory(ret)
            ret.atom_indices = self.atom_indices
        elif hasattr(ret, '_idx'):
            ret = ret.__subject__

        return ret

    def __reversed__(this):
        class ObjectIterator:
            def __init__(self):
                self.trajectory = this
                self.idx = len(this)
                self.length = 0

            def __iter__(self):
                return self

            def next(self):
                if self.idx > self.length:
                    self.idx -= 1
                    snapshot = self.trajectory[self.idx]
                    return snapshot.reversed
                else:
                    raise StopIteration()

        return ObjectIterator()

    def get_as_proxy(self, item):
        """
        Get an actual contained element

        This will also return lazy proxy objects and not the referenced ones
        as does __iter__, __reversed__ or __getitem__. Useful for faster access to the elements

        This is equal to use list.__getitem__(trajectory, item)

        Returns
        -------
        openpthsampling.Snapshot or openpathsampling.storage.objproxy.LoaderProxy
        """
        return list.__getitem__(self, item)


    def as_proxies(self):
        """
        Returns all contains all actual elements

        This will also return lazy proxy objects and not the references ones
        as does __iter__, __reversed__ or __getitme__. Useful for faster access to the elements

        Returns
        -------
        list of openpathsampling.Snapshot or openpathsampling.storage.objproxy.LoaderProxy
        """
        return list(self.iter_proxies())

    def iter_proxies(self):
        """
        Returns an iterator over all actual elements

        This will also return lazy proxy objects and not the references ones
        as does __iter__, __reversed__ or __getitme__. Useful for faster access to the elements

        Returns
        -------
        Iterator() over list of openpathsampling.Snapshot or openpathsampling.storage.objproxy.LoaderProxy

        """
        return list.__iter__(self)

    def __iter__(this):
        """
        Return an iterator over all snapshots in the storage

        This will always give real openpathsampling.Snapshot objects and never proxies to snapshots.
        If you prefer proxies (if available) use `.iteritems()`

        Parameters
        ----------
        iter_range : slice or None
            if this is not `None` it confines the iterator to objects specified
            in the slice

        Returns
        -------
        Iterator()
            The iterator that iterates the objects in the store

        """

        class ObjectIterator:
            def __init__(self):
                self.trajectory = this
                self.idx = 0
                self.length = len(this)

            def __iter__(self):
                return self

            def next(self):
                if self.idx < self.length:
                    obj = self.trajectory[self.idx]
                    self.idx += 1
                    return obj
                else:
                    raise StopIteration()

        return ObjectIterator()

    def __add__(self, other):
        t = Trajectory(self)
        t.extend(other)
        return t

    # =============================================================================================
    # PATH ENSEMBLE FUNCTIONS
    # =============================================================================================

    def summarize_by_volumes(self, label_dict):
        """Summarize trajectory based on number of continuous frames in volumes.

        This uses a dictionary of disjoint volumes: the volumes must be disjoint
        so that every frame can be mapped to one volume. If the frame maps to
        none of the given volumes, it returns the label None.

        Parameters
        ----------
        label_dict : dict
            dictionary with labels for keys and volumes for values

        Returns
        -------
        list of tuple
            format is (label, number_of_frames)
        """
        last_vol = None
        count = 0
        segment_labels = []
        # list.__iter__ for speed
        for frame in list.__iter__(self):
            in_state = []
            for key in label_dict.keys():
                vol = label_dict[key]
                if vol(frame):
                    in_state.append(key)
            if len(in_state) > 1:
                raise RuntimeError("Volumes given to summarize_by_volumes not disjoint")
            if len(in_state) == 0:
                current_vol = None
            else:
                current_vol = in_state[0]

            if last_vol == current_vol:
                count += 1
            else:
                if count > 0:
                    segment_labels.append( (last_vol, count) )
                last_vol = current_vol
                count = 1
        segment_labels.append( (last_vol, count) )
        return segment_labels

    def summarize_by_volumes_str(self, label_dict, delimiter="-"):
        """
        Return string version of the volumes visited by this trajectory.

        See `Trajectory.summarize_by_volumes` for details.

        Parameters
        ----------
        label_dict : dict
            dictionary with labels for keys and volumes for values
        delimiter : string (default "-")
            string used to separate volumes in output

        Returns
        -------
        string
            order in which this trajectory visits the volumes in
            `label_dict`, separated by the `delimiter`
        """
        summary = self.summarize_by_volumes(label_dict)
        return delimiter.join([str(s[0]) for s in summary])

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
            for snapshot_index in range(1, nsnapshots - 1):
                H += self[snapshot_index].kinetic_energy
        else:
            H = 0

        return H

    def computeActivity(self, atom_indices=None):
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
        K = 0.0

        if atom_indices is None:
            atom_indices = slice(None)

        for frame_index in range(nframes - 1):
            # Compute displacement of all atoms.
            delta_r = self[frame_index + 1].coordinates - self[frame_index].coordinates
            # Compute contribution to activity K.
            K += ((delta_r[atom_indices, :] / u.nanometers) ** 2).sum()

        return K * (u.nanometers ** 2)

    def logEquilibriumTrajectoryProbability(self):
        """
        Compute the (temperatureless!) log equilibrium probability

        Up to an unknown additive constant of an unbiased trajectory evolved
        according to Verlet dynamics with Andersen thermostatting.

        Parameters
        ----------
        trajectory : openpathsampling.Trajectory
            the trajectory

        Returns
        -------        
        log_q : float
            the log equilibrium probability of the trajectory divided by the
            inverse temperature beta
        
        NOTES
        -----
        This might be better places into the trajectory class. The trajectory
        should know the system and ensemble? and so it is not necessarily
        TPS specific

        """

        nsnapshots = len(self)
        log_q = - self[0].total_energy
        for snapshot_index in range(1, nsnapshots - 1):
            log_q += - self[snapshot_index].kinetic_energy

        return log_q

    # =============================================================================================
    # ANALYSIS FUNCTIONS
    # =============================================================================================

    def is_correlated(self, other):
        """
        Checks if two trajectories share a common snapshot

        Parameters
        ----------
        other : openpathsampling.Trajectory
            the second trajectory to check for common snapshots

        Returns
        -------
        bool
            returns True if at least one snapshot appears in both trajectories
        """

        # if hasattr(self, 'idx') and hasattr(other, 'idx'):
        #     shared_store = set(self.idx.keys()) & set(other.idx.keys())
        #     # both are saved so use the snapshot idx as identifiers
        #     if len(shared_store) > 0:
        #         storage = list(shared_store)[0]
        #         t1id = storage.trajectories.snapshot_indices(self.idx[storage])
        #         t2id = storage.trajectories.snapshot_indices(other.idx[storage])
        #         return bool(set(t1id) & set(t2id))

        # Use some fallback
        return bool(self.shared_configurations(other))

    def shared_configurations(self, other):
        """
        Returns a set of shared snapshots

        Parameters
        ----------
        other : Trajectory()
            the second trajectory to use

        Returns
        -------
        set of Snapshot()
            the set of common snapshots
        """
        return set([snap.configuration for snap in self]) & set([snap.configuration for snap in other])

    def shared_subtrajectory(self, other):
        """
        Returns a subtrajectory which only contains frames present in other

        Parameters
        ----------
        other : Trajectory()
            the second trajectory to use

        Returns
        -------
        Trajectory
            the shared subtrajectory
        """
        shared = self.shared_configurations(other)
        return Trajectory([snap for snap in self if snap.configuration in shared])

    # =============================================================================================
    # UTILITY FUNCTIONS
    # =============================================================================================

    def subset(self, atom_indices):
        """
        Reduce the view of the trajectory to a subset of atoms specified.

        This is only a view, no data will be changed or copied.
        
        Returns
        -------        
        trajectory : openpathsampling.Trajectory
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

        # TODO: To remove the dependency of the dynamics engine we need to get the information
        # TODO: about the solute_indices from somewhere else, preferrably the topology?

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

    def md(self, topology=None):
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

        velocities = \
            u.Quantity(np.zeros(mdtrajectory.xyz[0].shape, dtype=np.float32), u.nanometers / u.picoseconds)

        zero_momentum = paths.Momentum(velocities=velocities)

        for frame_num in range(len(mdtrajectory)):
            # mdtraj trajectories only have coordinates and box_vectors
            coord = u.Quantity(mdtrajectory.xyz[frame_num], u.nanometers)
            if mdtrajectory.unitcell_vectors is not None:
                box_v = u.Quantity(mdtrajectory.unitcell_vectors[frame_num],
                                   u.nanometers)
            else:
                box_v = None
            config = paths.Configuration(coordinates=coord, box_vectors=box_v)

            snap = paths.Snapshot(configuration=config, momentum=zero_momentum)
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
