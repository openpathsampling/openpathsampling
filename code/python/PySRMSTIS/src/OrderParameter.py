###############################################################
#| CLASS Order Parameter
###############################################################

from msmbuilder.metrics import RMSD
from trajectory import Trajectory
from snapshot import Snapshot
import numpy as np

class Cache(object):
    def __init__(self, fnc):
        '''
        Initializes an OrderParameter object that is essentially a function that maps a frame (Snapshot) within a trajectory (Trajectory) to a number. 
        
        Notes
        -----
        
        Does it have to be between 0 and 1?
        Do we want to save this in the netCDF?
        And if yes, do we want one big table for all OrderParameters or several ones. One for each. Second might be nicer.
        '''

        # Name needed to eventually store to netCDF. This is not yet done
        self._eval = fnc

        # We will use a cache since we assume the snapshots objects, the OrderParameter Instance and the ID# of snapshots to be immutable
        self.cache = []
        self.is_cached = set()
            
    def __getslice__(self, *args, **kwargs):
        """
        allows to use slicing and retains a Trajetory object!

        RETURNS
        
        trajectory (Trajectory) - the sliced trajectory
        
        NOTES
        
        This function is deprecated and will not be present in python 3 anymore
        
        """        
        sl = slice(*args, **kwargs)
        return self[range(sl.stop)[sl]]
                
    @property
    def size(self):
        return len(self.cache)
    
    def resize(self, size):
        if size > self.size:
            self.cache.extend([0.0 for i in range(self.size, size)])
    
    def in_cache(self, indices):
        return [idx for idx in indices if idx in self.is_cached]

    def not_in_cache(self, indices):
        return [idx for idx in indices if idx not in self.is_cached]   
            
    def __getitem__(self, index):
        """
        Adds the possibility to use indexing and retains a Trajetory object!
        
        RETURNS
        
        trajectory (Trajectory) - the indexed trajectory
                
        """        
                
        # Allow for numpy style of selecting several indices using a list as index parameter
        if type(index) is list:
            max_idx = max(index) + 1
            self.resize(max_idx)
            no_cache = self.not_in_cache(index)                        
        else:
            self.resize(index + 1)
            no_cache = self.not_in_cache([index])

        # Add not yet cached data            
        if len(no_cache) > 0:
            result = self._eval(no_cache)
        
            for no, res in enumerate(result):
                self.cache[no_cache[no]] = res

            self.is_cached.update(no_cache)

        if type(index) is list:
            ret = [self.cache[idx] if idx > 0 else None for idx in index]                        
        else:                        
            ret = self.cache[index]  
                
        return ret
    
class SnapshotCache(Cache):
    
    def __init__(self, name, fnc):
        super(SnapshotCache, self).__init__(fnc = fnc)
        
        self.name = name

    def fill(self):
        '''
        Compute all distances for all snapshots to be used later using the cache. 
        
        Notes
        -----
        
        Pick another name?
        Make sure that all snapshots are saved. Otherwise we cannot cache them!
        
        '''

        print self[1:Snapshot.load_number() + 1]


    def _init_netcdf(self, storage):
        '''
        initializes the associated storage to save a specific order parameter in it
        '''           
        # save associated storage in class variable for all Snapshot instances to access
        
#        ncgrp = storage.ncfile.createGroup('snapshot')
        
        self.storage = storage
        ncgrp = storage.ncfile
        
        # define dimensions used in snapshots
#        size_dimension = 'op_size_' + self.name
#        ncgrp.createDimension(size_dimension, self.size)                       # unlimited number of snapshots
        
        # define variables for OrderParameters
        var_name = 'order_parameter_' + self.name
        ncvar_order_parameter          = ncgrp.createVariable(var_name, 'f', ('snapshot', 'scalar'))

        # Define units for snapshot variables.
        setattr(ncvar_order_parameter, 'units', 'None')
        
        # Define long (human-readable) names for variables.
        setattr(ncvar_order_parameter,    "long_name", "orderparameter[snapshot][index] is the orderparameter of snapshot 'snapshot'.")
        
    def save(self):
        """
        Save the current state of the cache to the storage
        """
        
        ncfile = self.storage.ncfile

        # Make sure snapshots are stored and have an index and then add the snapshot index to the trajectory
        var_name = 'order_parameter_' + self.name
        indices = list(self.is_cached)
        ncfile.variables[var_name][:,0] = np.array(self.cache)
        ncfile.variables[var_name + '_idx'] = np.array(indices) 
        
        
    def load(self):
        '''
        Restores the cache from the storage        
        '''
        var_name = 'order_parameter_' + self.name
        self.cache = self.storage.ncfile.variables[var_name][:,0].astype(np.float).copy().as_list()
        self.is_cached = set(self.storage.ncfile.variables[var_name + '_idx'].astype(np.int).copy().as_list())  

class OrderParameter(object):
    def __init__(self, name):
        '''
        Initializes an OrderParameter object that is essentially a function that maps a frame (Snapshot) within a trajectory (Trajectory) to a number. 
        
        Notes
        -----
        
        Does it have to be between 0 and 1?
        Do we want to save this in the netCDF?
        And if yes, do we want one big table for all OrderParameters or several ones. One for each. Second might be nicer.
        '''

        # Name needed to eventually store to netCDF. This is not yet done
        self.name = name

        # We will use a cache since we assume the snapshots objects, the OrderParameter Instance and the ID# of snapshots to be immutable
        self.use_cache = True
        self.use_storage = False
        
        # Run the OrderParameter on the initial snapshot to get the size
#        self.size = len(self._eval(0))
        self.size = 1
        
        if self.use_cache:
            self.cache = SnapshotCache(name = self.name, fnc = self._eval)
        else:
            self.cache = None
            
    def _eval(self, indices, path = None):
        '''
        Actual evaluation of indices of snapshots. This needs to be implemented.
        '''
        pass

    def __call__(self, snapshot):
        '''
        Calculates the actual order parameter by making the OrderParameter object into a function
        
        Parameters
        ----------

        snapshot (Snapshot)
            snapshot object used to compute the lambda value
        
        Notes
        -----
        
        '''
        return self._assign(snapshot)        
    

    
    def _assign(self, snapshots):
        '''
        Assign a single snapshot or a trajectory.
        
        Parameters
        ----------
        
        snapshot : Snapshot
            A list of or a single snapshot that should be assigned
            
        Notes
        -----
        This calls self._eval with the relevant snapshots.
        '''
        single = False
        
        if type(snapshots) is Snapshot:
            traj = Trajectory([snapshots])
            single = True
        else:
            traj = snapshots

        traj_indices = traj.indices()
        
        if self.use_cache:
            # pick all non-stored snapshots and compute these anyway without storage
            no_index = [ i for i,idx in enumerate(traj_indices) if idx == 0]
            
            d = self.cache[traj_indices]

            if len(no_index):
                # compute ones that cannot be cached
                d_no_index = self._eval(no_index)             
    
                # add unknowns
                for ind, dist in enumerate(d_no_index):
                    d[no_index[ind]] = dist
                                                
        else: 
            d = self._eval(traj_indices)

        if single:
            return d[0]
        else:
            return d    
    
    @staticmethod
    def _scale_fnc(mi, ma):
        '''
        Helper function that returns a function that scale values in a specified range to a range between 0 and 1.
        
        Parameters
        ----------
        
        mi : float
            Minimal value. Corresponds to zero.
        ma : float
            Maximal value. Corresponds to one.
            
        Returns
        -------
        scale : function
            The scaling function of float -> float
        
        '''
        
        def scale(x):
            if x < mi:
                return 0.0
            elif x > ma:
                return 1.0
            else:
                return (x-mi) / (ma-mi)
            
        return scale
        
        
class OP_RMSD_To_Lambda(OrderParameter):
    def __init__(self, name, center, lambda_min, max_lambda, atom_indices = None):
        '''
        A OrderParameter that transforms the RMSD to a specific center to a lambda value between zero and one.
        
        Parameters
        ----------
        
        center : snapshot
            a trajectory snapshot that is used as the point to compute the RMSD to
        lambda_min : float
            rmsd value that corresponds to lambda zero
        max_lambda : float
            rmsd value that corresponds to lambda one
        atom_indices : list of integers (optional)
            a list of integers that is used in the rmsd computation. Usually solvent should be excluded
        
        '''
        
        super(OP_RMSD_To_Lambda, self).__init__(name)
                
        self.atom_indices = atom_indices
        self.center = center
        self.min_lambda = lambda_min
        self.max_lambda = max_lambda
        
        self.size = 1
                
        # Generate RMSD metric using only the needed indices. To save memory we crop the read snapshots from the database and do not use a cropping RMSD on the full snapshots
        self.metric = RMSD(None)
        self.generator = self.metric.prepare_trajectory(Trajectory([center]).subset(self.atom_indices).md())  
        return 
    
    ################################################################################
    ##  Actual computation of closest point using RMSD
    ################################################################################    
    
    def _eval(self, indices):
        ptraj = self.metric.prepare_trajectory(Trajectory.from_indices(indices).subset(self.atom_indices).md())
        results = self.metric.one_to_all(self.generator, ptraj, 0)
                
        return map(self._scale_fnc(self.min_lambda, self.max_lambda), results )
    
if __name__ == '__main__':
    def ident(indices):
        if type(indices) is list:
            return [float(i) for i in indices]
        else:
            return float(indices)
        
    s = SnapshotCache(name = 'TestList', fnc = ident)
    
    print s[10]
    print s[5]
    print s[5]
    print s.cache
    print s[1:10]
    print s.cache