#!/usr/bin/env python
from __future__ import print_function, absolute_import, division
import numpy as np
from msmbuilder import clustering
from msmbuilder import metrics
import logging

logger = logging.getLogger('msmbuilder.scripts.Cluster')

################################################################################

def assign_in_memory(metric, generators, storage, atom_indices=None):

    traj = storage.all_snapshot_coordinates_as_array(atom_indices)
    
    n_frames = len(traj)

    assignments = -1 * np.ones(n_frames, dtype='int')
    distances = -1 * np.ones(n_frames, dtype='float32')

    pgens = metric.prepare_trajectory(generators)    
    ptraj = metric.prepare_trajectory(traj)

    for j in xrange(len(traj)):
        d = metric.one_to_all(ptraj, pgens, j)
        assignments[j] = np.argmin(d)
        distances[j] = d[assignments[j]]

    return assignments, distances

    def update_cluster_from_storage(metric):
                
        if isinstance(metric, metrics.RMSD):
            atom_indices = metric.atomindices
            metric.atomindices = None # probably bad...
            logger.info('RMSD metric - loading only the atom indices required')
        else:
            atom_indices = None
    
        
        traj = storage.all_snapshot_coordinates_as_array(atom_indices)
         
        ptrajs = None
        
        
        
        
        clusterer = clustering.HybridKMedoids(metric, trajectories=[traj],
            prep_trajectories=ptrajs, k=args.hybrid_num_clusters,
            distance_cutoff=args.hybrid_distance_cutoff,
            local_num_iters=args.hybrid_local_num_iters,
            global_num_iters=args.hybrid_global_iters,
            too_close_cutoff=args.hybrid_too_close_cutoff,
            ignore_max_objective=args.hybrid_ignore_max_objective)
    
        gen_inds = clusterer.get_generator_indices()
        
        
def get_cluster(metric, generators, state, atom_indices=None):

    pgens = metric.prepare_trajectory(generators)    
    ptraj = metric.prepare_trajectory([state])

    d = metric.one_to_all(ptraj, pgens, 0)
    assignment = np.argmin(d)
    distance = d[assignment]

    return assignment, distance
