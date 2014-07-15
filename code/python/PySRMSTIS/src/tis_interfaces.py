'''
Created on 07.07.2014

@author: JH Prinz

NOTES

Main file copied from Cluster.py from MSMBuilder
'''

from __future__ import print_function, absolute_import, division
import numpy as np

from pymbar import MBAR

from msmbuilder.metrics import (RMSD)

class Object(object):
            pass

from trajectory import Trajectory

class TISInterfaces(object):
    '''
    Hold the TIS Cores and Inferfaces description and can decide where a snapshot is located, etc...
    
    NOTES
    
    Is the core / interface numbering unique? I don't think so, but then we need to add this
    '''

    def __init__(self, storage, atom_indices, generators, distances):
        '''
        Constructor
        '''
        
        self.storage = storage
        self.generators = generators
        self.lambdas = distances
        self.metric = None
        self.atom_indices = atom_indices
    
    ################################################################################

    def assign_from_indices(self, indices):
        return (self.cores[indices], self.interfaces[indices], self.distances[indices])
        
    def assign_all_trajectories(self):
        return [ self.assign_from_indices(indices) for indices in self.storage.all_trajectory_indices() ]
    
    def assign_storage(self):
        
        traj = self.storage.all_snapshot_coordinates_as_mdtraj( self.atom_indices )
        
        n_frames = len(traj)
    
        assignments = -1 * np.ones(n_frames, dtype='int')
        interfaces = -1 * np.ones(n_frames, dtype='int')
        distances = -1 * np.ones(n_frames, dtype='float32')
        
        pgens = self.metric.prepare_trajectory(self.generators)    
        ptraj = self.metric.prepare_trajectory(traj)

        for j in xrange(len(traj)):
            d = self.metric.one_to_all(ptraj, pgens, j)
            assignments[j] = np.argmin(d)
            distances[j] = d[assignments[j]]
            l_lambdas = self.lambdas[assignments[j]]
            
            if distances[j] < l_lambdas[-1]:
                interfaces[j] = len([l for l in l_lambdas if l < distances[j]]) 
            else:
                assignments[j] = -1
                
        self.cores = assignments
        self.interfaces = interfaces
        self.distances = distances

        return

    def assign_trajectory(self, traj):
        indices = traj.indices()
        return self.assign_from_indices(indices)
            
    def assign_snapshot(self, snapshot):
        
        traj = Trajectory([snapshot]).subset(self.atom_indices).md()
            
        assignments = -1 
        interfaces = -1
        distances = -1
        
        pgens = self.metric.prepare_trajectory(self.generators)    
        ptraj = self.metric.prepare_trajectory(traj)

        d = self.metric.one_to_all(ptraj, pgens, 0)

        assignments = np.argmin(d)
        distances = d[assignments]
        l_lambdas = self.lambdas[assignments]
        
        if distances < l_lambdas[-1]:
            interfaces = len([l for l in l_lambdas if l < distances]) 
        else:
            assignments = -1
                
        return assignments, interfaces, distances
    
    def assign(self, traj, atom_indices=None):

        n_frames = len(traj)
    
        assignments = -1 * np.ones(n_frames, dtype='int')
        interfaces = -1 * np.ones(n_frames, dtype='int')
        distances = -1 * np.ones(n_frames, dtype='float32')
    
        pgens = self.metric.prepare_trajectory(self.generators)    
        ptraj = self.metric.prepare_trajectory(traj.md( atom_indices ))

        for j in xrange(len(traj)):
            d = self.metric.one_to_all(ptraj, pgens, j)
            assignments[j] = np.argmin(d)
            distances[j] = d[assignments[j]]
            interfaces[j] = np.argmax(self.distances[self.distances < distances[j]])

        return assignments, interfaces, distances

    def in_interface(self, snapshot):
        c, i, d = self.assign_snapshot(snapshot)
        return c != -1 and i != -1       # hit a core
    
    def within_interface(self, idx):
        def stopper(snapshot):
            c, i, d = self.assign_snapshot(snapshot)
            return c != -1 and i == idx      # hit a core
        
        return stopper

    def in_core(self, snapshot):
        c, i, d = self.assign_snapshot(snapshot)
        return c != -1 and i == 0      # hit a core
    
    def create_rmsd_metric(self, atom_indices = None):
        self.metric = RMSD(atom_indices)  # , omp_parallel=args.rmsd_omp_parallel)
                
    def split_into_connections(self, traj):
        c, i, d = self.assign_trajectory(traj)
        
        l_traj = []
        traj_mode = 0
        t = Trajectory()
        first = None
        for index, snapshot in enumerate(traj):
            if (c[index] != -1 and i[index]==0 and traj_mode == 0):
                # core set
                first = snapshot
            elif (i[index] != 0  and first is not None):
                # in void
                t.append(snapshot)
                traj_mode = 1                
            elif (c[index] != -1 and i[index]==0 and traj_mode == 1):
                t.insert(0, first)
                t.append(snapshot)
                l_traj.append(t)
                t = Trajectory()
                first = snapshot
                traj_mode = 0
        return l_traj
    
    def split_into_connections_indices(self, indices):
        c, i, d = self.assign_from_indices(indices)
        
        l_traj = []
        t = []
        traj_mode = 0
        first = None
        for index in range(len(indices)):
            if (c[index] != -1 and i[index]==0 and traj_mode == 0):
                # core set
                first = index
            elif (i[index] != 0  and first is not None):
                # in void
                t.append(indices[index])
                traj_mode = 1                
            elif (c[index] != -1 and i[index]==0 and traj_mode == 1):
                t.insert(0, first)
                t.append(indices[index])
                l_traj.append(t)
                t = Trajectory()
                first = indices[index]
                traj_mode = 0
        return l_traj
    
    def _get_TIS_Info_from_connection(self, indices):
        '''
        Computes the first core, the last core and a list of crossed interfaces from the first state as well as the maximum interface. 
        
        NOTES:
        
        The list of crossed interfaces given a maximal interface - does it need to contain all interfaces that are smaller then the max? And especially, if it goes from I -> J
        does it then (unless cores are very close) need to cross all interfaces?
        '''
        c, i, d = self.assign_from_indices(indices)
                
        if c[0] != -1 and i[0] == 0:
            first = c[0]
        else:
            # ERROR: Does not start in core
            first = -1    

        if c[-1] != -1 and i[-1] == 0:
            last = c[-1]
        else:
            print(c, i, d)
            # Error: Does not end in core
            last = -1    

        
        n_interfaces = self.lambdas[first].shape[0]
        i_matrix = np.zeros(n_interfaces, dtype='int')
                    
        for idx in range(1,len(c) - 1):
            if c[idx] == first and i[idx] != -1:
                i_matrix[i[idx]] = 1
        
#        print(first, last)
#        print(i_matrix)
        max_lambda_I = np.max(np.array(range(n_interfaces))[(i_matrix == 1)])
#        print(max_lambda_I)
                
        return first, last, max_lambda_I, (i_matrix == 1)

    def compute_N_IJ_l_kI(self):
        n_cores = len(self.generators)
        n_interfaces = self.lambdas.shape[1] #for now assume that all cores have the same number of interfaces
        n_matrix = np.zeros((n_cores, n_cores, n_interfaces, n_interfaces), dtype='int')
        
        all_traj_indices = self.storage.all_trajectory_indices()
        
        for indices in all_traj_indices:
            connectors = self.split_into_connections_indices(indices)
            for connection in connectors:
                f, l, m, k = self._get_TIS_Info_from_connection(connection)
                n_matrix[f,l,k,m] += 1
                    
        return n_matrix