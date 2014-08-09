###############################################################
#| HEADER
###############################################################

from msmbuilder.metrics import RMSD
from trajectory import Trajectory
from snapshot import Snapshot

class OrderParameter(object):
    def __init__(self):
        '''
        Initializes an OrderParameter object that is essentially a function that maps a frame (Snapshot) within a trajectory (Trajectory) to a number in [0,1] 
        
        NOTES
        
        Does it have to be between 0 and 1?
        '''
        pass
    
    

    def __call__(self, frame, path = None):
        '''
        Calculates the actual order parameter by making the OrderParameter object into a function
        
        PARAMETERS
        frame (Snapshot) - snapshot object used to compute the lambda value
        
        OPTIONAL PARAMETERS
        path (Trajectory) - trajectory object that might be used in the Order Parameter
        
        NOTES
        
        '''
        
        
        pass

class OP_RMSD_To_Snaphost(OrderParameter):
    def __init__(self, center, min_lambda, max_lambda, atom_indices = None):
        '''
        
        '''
                
        self.atom_indices = atom_indices
        self.center = center
        self.min_lambda = min_lambda
        self.max_lambda = max_lambda
        
        # We will use a cache since we assume the snapshots objects, the OrderParameter Instance and the ID# of snapshots to be immutable
        self.cache = []
        
        # Generate RMSD metric using only the needed indices. To save memory we crop the read snapshots from the database and not use a cropping RMSD on the full snapshots
        self.metric = RMSD(None)
        self.generator = self.metric.prepare_trajectory(Trajectory([center]).md())  
        return 
    
    ################################################################################

    
    def assign_all(self):
        '''
        Compute all distances for all snapshots to be used later. 
        
        NOTES
        
        Pick another name?
        '''
        pass
    
    
    def _assign(self, snapshots):
        '''
        Assign a single snapshot or a trajectory
        '''
        
        use_cache = True
        
        single = False
        
        if type(snapshots) is Snapshot:
            traj = Trajectory([snapshots])
            single = True
        else:
            traj = snapshots
        
        if use_cache:
            # Remove precalculated distances
            
            traj_indices = Trajectory().indices()
            max_cached = len(self.cache)

            # pick all snapshots that have a number, others cannot be stored
            with_index = list(set([ idx for i,idx in enumerate(traj_indices) if idx > 0]))
            
            # pick all non-stored snapshots and compute these anyway without storage
            no_index = [ idx for i,idx in enumerate(traj_indices) if idx == 0]
            
            # limit to ones that could be cached
            in_cache = [ idx for i in with_index if max_cached >= idx]
            
            # check if it has actually been stored
            in_cache = [ i for i in in_cache if self.cache[traj_indices[i]] >= 0.0 ]
            
            # idx that should be added to the cache. Note that only snapshots with an ID can be added to the cache
            to_cache = [ i for i in with_index if i not in in_cache ]
                        
            ptraj = self.metric.prepare_trajectory(Trajectory.load_indices(to_cache).subset(self.atom_indices).md())
            d_to_cache = self.metric.one_to_all(self.generator, ptraj, 0)
            
            # add to cache
            for ind, dist in enumerate(d_to_cache):
                self.cache[to_cache[ind]] = dist

            ptraj = self.metric.prepare_trajectory(Trajectory.load_indices(no_index).subset(self.atom_indices).md())
            d_no_index = self.metric.one_to_all(self.generator, ptraj, 0)

            d = []

            # use cache to get distances
            for i, idx in enumerate(traj_indices):
                if idx > 0:
                    d[i] = self.cache[idx]

            # add unknowns
            for ind, dist in enumerate(d_no_index):
                d[no_index[ind]] = dist
                                                
        else: 
            ptraj = self.metric.prepare_trajectory(traj.subset(self.atom_indices).md())
            d = self.metric.one_to_all(self.generator, ptraj, 0)

        if single:
            return d[0]
        else:
            return d


    def __call__(self, snapshot, path = None):
        '''
        Override function of the OrderParameter and return the actual value
        '''
        
        return self._assign(snapshot)