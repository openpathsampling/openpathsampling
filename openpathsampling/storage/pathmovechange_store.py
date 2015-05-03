from object_storage import ObjectStore
from openpathsampling.pathmovechange import PathMoveChange

from object_storage import ObjectStore
from openpathsampling.pathmovechange import PathMoveChange
from openpathsampling.storage.object_storage import func_update_object

import time

class PathMoveChangeStore(ObjectStore):
    def __init__(self, storage):
        super(PathMoveChangeStore, self).__init__(storage, PathMoveChange, json=False, load_partial=True)

        self.set_variable_partial_loading('details', self.update_details)
        self.set_variable_partial_loading('parent', self.update_parent)
        self.set_variable_partial_loading('mover', self.update_mover)
        self.set_variable_partial_loading('ensemble', self.update_ensemble)


    def load_empty(self, idx):
        trajectory_idx = int(self.storage.variables['pathmovechange_trajectory_idx'][idx])
        replica_idx = int(self.storage.variables['pathmovechange_replica'][idx])
        valid=self.load_variable('pathmovechange_valid', idx)
        accepted=bool(self.load_variable('pathmovechange_accepted', idx))

        obj = PathMoveChange(
            trajectory=self.storage.trajectories[trajectory_idx],
            replica=replica_idx,
            valid=valid,
            accepted=accepted,
        )

        del obj.details
        del obj.ensemble
        del obj.mover
        del obj.parent

        return obj

    update_details = func_update_object('pathmovechange', 'details', '_details')
    update_parent = func_update_object('pathmovechange', 'parent', 'trajectories')
    update_mover = func_update_object('pathmovechange', 'mover', 'pathmovers')
    update_ensemble = func_update_object('pathmovechange', 'ensemble', 'ensembles')

    def save(self, pathmovechange, idx=None):
        if idx is not None:
            self.storage.trajectories.save(pathmovechange.trajectory)
            self.save_object('pathmovechange_trajectory', idx, pathmovechange.trajectory)

            self.storage.ensembles.save(pathmovechange.ensemble)
            self.save_object('pathmovechange_ensemble', idx, pathmovechange.ensemble)

            self.save_variable('pathmovechange_replica', idx, pathmovechange.replica)
            self.save_object('pathmovechange_parent', idx, pathmovechange.parent)
            self.save_object('pathmovechange_details', idx, pathmovechange.details)
            self.save_variable('pathmovechange_valid', idx, pathmovechange.valid)
            self.save_variable('pathmovechange_accepted', idx, pathmovechange.accepted)
            self.save_object('pathmovechange_pathmover', idx, pathmovechange.mover)

    def load(self, idx):
        '''
        Return a pathmovechange from the storage

        Parameters
        ----------
        idx : int
            index of the pathmovechange

        Returns
        -------
        pathmovechange : PathMoveChange
            the pathmovechange
        '''
        trajectory_idx = int(self.storage.variables['pathmovechange_trajectory_idx'][idx])
        ensemble_idx = int(self.storage.variables['pathmovechange_ensemble_idx'][idx])
        replica_idx = int(self.storage.variables['pathmovechange_replica'][idx])
        parent_idx = int(self.storage.variables['pathmovechange_parent_idx'][idx])
        details_idx = int(self.storage.variables['pathmovechange_details_idx'][idx])
        pathmover_idx = int(self.storage.variables['pathmovechange_pathmover_idx'][idx])
        valid=self.load_variable('pathmovechange_valid', idx)
        accepted=bool(self.load_variable('pathmovechange_accepted', idx))


        obj = PathMoveChange(
            trajectory=self.storage.trajectories[trajectory_idx],
            replica=replica_idx,
            ensemble=self.storage.ensembles[ensemble_idx],
            valid=valid,
            parent=self.storage.pathmovechanges[parent_idx],
            details=self.storage._details[details_idx],
            accepted=accepted,
            mover=self.storage.pathmovers[pathmover_idx]
        )

        return obj

    def by_ensemble(self, ensemble):
        return [ pathmovechange for pathmovechange in self.iterator() if pathmovechange.ensemble == ensemble ]

    def _init(self, units=None):
        super(PathMoveChangeStore, self)._init(units)

        # New short-hand definition
        self.init_variable('pathmovechange_trials_idx', 'index', 'sampleset',
            description="sampleset[sampleset][frame] is the sample index (0..nspanshots-1) of frame 'frame' of sampleset 'sampleset'.",
            variable_length = True,
            chunksizes=(1024, )
        )
        self.init_variable('pathmovechange_trials', 'index', chunksizes=(1, ))
        self.init_variable('pathmovechange_ensemble_idx', 'index', chunksizes=(1, ))
        self.init_variable('pathmovechange_replica', 'index', chunksizes=(1, ))
        self.init_variable('pathmovechange_parent_idx', 'index', chunksizes=(1, ))
        self.init_variable('pathmovechange_valid', 'index', chunksizes=(1, ))
        self.init_variable('pathmovechange_details_idx', 'index', chunksizes=(1, ))
        self.init_variable('pathmovechange_accepted', 'bool', chunksizes=(1, ))
        self.init_variable('pathmovechange_pathmover_idx', 'index', chunksizes=(1, ))



class PathMoveChangeStore(ObjectStore):
    def __init__(self, storage):
        super(PathMoveChangeStore, self).__init__(storage, PathMoveChange, has_uid=False, nestable=True, load_partial=True)
        self.set_variable_partial_loading('details', self.update_details)

    def load_empty(self, idx):
        trajectory_idx = int(self.storage.variables['pathmovechange_trajectory_idx'][idx])
        ensemble_idx = int(self.storage.variables['pathmovechange_ensemble_idx'][idx])
        replica_idx = int(self.storage.variables['pathmovechange_replica'][idx])
        step=self.load_variable('pathmovechange_step', idx)

        if step < 0:
            step = None


        obj = PathMoveChange(
            trajectory=self.storage.trajectories.load(trajectory_idx),
            replica=replica_idx,
            ensemble=self.storage.ensembles.load(ensemble_idx),
            step=step
        )

        del obj.details

        return obj

    @staticmethod
    def update_details(obj):
        storage = obj._origin

        idx = obj.idx[storage]
        details_idx = int(storage.variables['pathmovechange_details_idx'][idx])
        details = storage.movedetails.load(details_idx)

        obj.details = details

    def load(self, idx):
        '''
        Return a pathmovechange from the storage

        Parameters
        ----------
        idx : int
            index of the pathmovechange

        Returns
        -------
        pathmovechange : PathMoveChange
            the pathmovechange
        '''
        trajectory_idx = int(self.storage.variables['pathmovechange_trajectory_idx'][idx])
        ensemble_idx = int(self.storage.variables['pathmovechange_ensemble_idx'][idx])
        replica_idx = int(self.storage.variables['pathmovechange_replica'][idx])
        details_idx = int(self.storage.variables['pathmovechange_details_idx'][idx])
        step=self.load_variable('pathmovechange_step', idx)


        obj = PathMoveChange(
            trajectory=self.storage.trajectories.load(trajectory_idx),
            replica=replica_idx,
            ensemble=self.storage.ensembles.load(ensemble_idx),
            details=self.storage.movedetails.load(details_idx),
            step=step
        )

        return obj

    def by_ensemble(self, ensemble):
        return [ pathmovechange for pathmovechange in self.iterator() if pathmovechange.ensemble == ensemble ]

    def _init(self):
        super(PathMoveChangeStore, self)._init()

        # New short-hand definition
        self.init_variable('pathmovechange_trajectory_idx', 'index', chunksizes=(1, ))
        self.init_variable('pathmovechange_ensemble_idx', 'index', chunksizes=(1, ))
        self.init_variable('pathmovechange_replica', 'index', chunksizes=(1, ))
        self.init_variable('pathmovechange_details_idx', 'index', chunksizes=(1, ))
        self.init_variable('pathmovechange_step', 'index', chunksizes=(1, ))

class PathMoveChangeSetStore(ObjectStore):

    def __init__(self, storage):
        super(PathMoveChangeSetStore, self).__init__(storage, PathMoveChangeSet, json=False)

    def save(self, pathmovechange_set, idx=None):
        # Check if all pathmovechanges are saved
        map(self.storage.pathmovechanges.save, pathmovechange_set)

        values = self.list_to_numpy(pathmovechange_set, 'pathmovechange')
        self.storage.variables['pathmovechangeset_pathmovechange_idx'][idx] = values

        self.storage.pathmovechanges.save(pathmovechange_set.movepath)
        self.save_object('pathmovechangeset_movepath', idx, pathmovechange_set.movepath)


    def pathmovechange_indices(self, idx):
        '''
        Load pathmovechange indices for pathmovechange_set with ID 'idx' from the storage

        Parameters
        ----------
        idx : int
            ID of the pathmovechange_set

        Returns
        -------
        list of int
            list of pathmovechange indices
        '''

        # get the values
        values = self.storage.variables['pathmovechangeset_pathmovechange_idx'][idx]

        # typecast to integer
        return self.list_from_numpy(values, 'index')

    def load(self, idx):
        '''
        Return a pathmovechange_set from the storage

        Parameters
        ----------
        idx : int
            index of the pathmovechange_set (counts from 0)

        Returns
        -------
        pathmovechange_set
            the pathmovechange_set
        '''

        movepath_idx = int(self.storage.variables['pathmovechangeset_movepath_idx'][idx])

        values = self.storage.variables['pathmovechangeset_pathmovechange_idx'][idx]

        # typecast to pathmovechange
        pathmovechanges = self.list_from_numpy(values, 'pathmovechanges')
        pathmovechange_set = PathMoveChangeSet(pathmovechanges, movepath=self.storage.pathmovechanges.load(movepath_idx))

        return pathmovechange_set

    def _init(self, units=None):
        """
        Initialize the associated storage to allow for pathmovechangeset storage

        """
        super(PathMoveChangeSetStore, self)._init()

        self.init_variable('pathmovechangeset_pathmovechange_idx', 'index', 'pathmovechangeset',
            description="pathmovechangeset[pathmovechangeset][frame] is the pathmovechange index (0..nspanshots-1) of frame 'frame' of pathmovechangeset 'pathmovechangeset'.",
            variable_length = True,
            chunksizes=(1024, )
        )

        self.init_variable('pathmovechangeset_movepath_idx', 'index', chunksizes=(1, ))
