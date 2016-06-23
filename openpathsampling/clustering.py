'''
Created on 07.07.2014

@author: JH Prinz

NOTES

Main file copied from Cluster.py from MSMBuilder
'''

from __future__ import print_function, absolute_import, division
import numpy as np
from msmbuilder.cluster import clustering

from msmbuilder.metrics import RMSD

from openpathsampling import Trajectory

class Object(object):
            pass

class VoronoiTesselation(object):
    '''
    Hold the MSM Clustering description and the associated state assignments to decide if a trajectory has hit a core
    
    Notes
    -----
    
    
    '''

    def __init__(self):
        '''
        Create an empty Voronoi Tessalation object
        
        NOTES
        
        '''
        
        self.storage = None
        self._generator = None
        self.metric = RMSD(None)
        self.atom_indices = None
        self.snapshot_distances = None
        self.snapshot_indices = None
    
    ################################################################################

    def update_cluster_from_storage(self):
        '''
        Update the set of generators from the associates trajectory storage
        
        Notes
        -----
        
        '''
        traj = self.storage.all_snapshot_coordinates_as_mdtraj( self.atom_indices )
        
        args = Object()
        
        args.hybrid_local_num_iters = 50
        args.hybrid_global_iters = 0
        args.hybrid_ignore_max_objective = False
        args.hybrid_too_close_cutoff = 0.0001
        args.hybrid_num_clusters = self.n_centers
        args.hybrid_distance_cutoff = None

        ptrajs = None
    
        clusterer = clustering.HybridKMedoids(
            self.metric, 
            trajectories=traj,
            prep_trajectories=ptrajs, 
            k=args.hybrid_num_clusters,
            distance_cutoff=args.hybrid_distance_cutoff,
            local_num_iters=args.hybrid_local_num_iters,
            global_num_iters=args.hybrid_global_iters,
            too_close_cutoff=args.hybrid_too_close_cutoff,
            ignore_max_objective=args.hybrid_ignore_max_objective
        )
    
        gen_inds = clusterer.get_generator_indices()

        self.generators_indices = gen_inds
        self.generators = traj[gen_inds]
        
        return
    
    @property
    def size(self):
        '''
        Return the number of generators used in the tesselation
        
        Returns
        -------
        length : int
            number of generators

        '''

        return len(self.generators)

    def assign_storage(self):
        '''
        Assign all snapshots in the associates trajectory storage to the generators
                
        Notes
        -----
        This allows later to access everything fast
        '''
        
        traj = self.storage.all_snapshot_coordinates_as_mdtraj( self.atom_indices )
        
        n_frames = len(traj)
    
        assignments = -1 * np.ones(n_frames, dtype='int')
        distances = -1 * np.ones(n_frames, dtype='float32')
    
        pgens = self.metric.prepare_trajectory(self.generators)    
        ptraj = self.metric.prepare_trajectory(traj)

        for j in xrange(len(traj)):
            d = self.metric.one_to_all(ptraj, pgens, j)
            assignments[j] = np.argmin(d)
            distances[j] = d[assignments[j]]
        
        self.snapshot_indices = np.array(assignments, dtype='int')
        self.snapshot_distances = np.array(distances, dtype='float')
        
        return
    
    def assign_all_trajectories(self):
        '''
        Assign all trajectories in the associates trajectory storage to the generators
        
        Returns
        -------
        clusterlist : list of int
            list of cluster IDs
        
        Notes
        -----        
        This needs assign_storage() to be run before!
        
        '''

        return [ self.snapshot_indices[t] for t in self.storage.all_trajectory_indices() ]

    def assign_index_trajectory(self, indices):
        '''
        Assign snapshots with IDs indices to the generators
        
        Returns
        -------
        clusterlist (list of int) - list of cluster IDs
        
        NOTES
        
        '''

        if self.snapshot_indices is not None:
            return self.snapshot_indices[indices]
        else:
            return None
            
    def assign_snapshot(self, snapshot):
        '''
        Assign a single snapshot to the cluster centers
        
        Returns
        -------
        assignment : int
            cluster IDs
        distance: float
            distance to cluster _generator in measure of the metric (RMSD)
        
        '''
        
        assignments, distances = self.assign(Trajectory([snapshot]))    
        return assignments[0], distances[0]

    def assign(self, traj, recalc = False):
        '''
        Assign a Trajectory object to the cluster
        
        Parameters
        ----------        
        traj : Trajectory
            trajectory to be clustered
        recalc : bool
            forces a calculation of the cluster center and not using the cached assignments (Default False)
        
        RETURNS
        
        assignments (numpy.array(n_frames, dtype='int')) - array of cluster IDs
        distances (numpy.array(n_frames, dtype='float')) - distances to cluster _generator in measure of the metric (RMSD)
                
        '''

        n_frames = len(traj)
            
        if self.snapshot_indices is not None:
            # We do not check if the snapshot_indices are properly updated!
            # Checking might be too expensive
            indices = traj.indices()
            assignments = self.snapshot_indices[indices]
            distances = self.snapshot_distances[indices]
            
        else:
            assignments = -1 * np.ones(n_frames, dtype='int')
            distances = -1 * np.ones(n_frames, dtype='float32')

            pgens = self.metric.prepare_trajectory( self.generators )    
            ptraj = self.metric.prepare_trajectory( traj.subset(self.atom_indices).md() )
    
            for j in xrange(len(traj)):
                d = self.metric.one_to_all(ptraj, pgens, j)
                assignments[j] = np.argmin(d)
                distances[j] = d[assignments[j]]

        return assignments, distances