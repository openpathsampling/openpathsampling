from object_storage import ObjectStorage

class Origin(object):
    """
    A Move is the return object from a PathMover and contains all information about the move, initial trajectories,
    new trajectories (both as references). Might move several trajectories at a time (swapping)

    Notes
    -----
    Should contain inputs/outputs and success (accepted/rejected) as well as probability to succeed.
    """

    cls = 'origin'

    def __init__(self, name=None, final=None, result=None, inputs=None, mover=None, acceptance=None, success=None, options=None, ensemble=None):
        self.name = name
        self.inputs = inputs
        self.mover = mover
        self.final = final
        self.acceptance = acceptance
        self.success = success
        self.options = options
        self.ensemble = ensemble
        self.result = result

        self.idx = dict()

class OriginStorage(ObjectStorage):
    def __init__(self, storage):
        super(OriginStorage, self).__init__(storage, Origin)

    def save(self, origin, idx=None):
        """
        Add the current state of the origin in the database. If nothing has changed then the origin gets stored using the same snapshots as before. Saving lots of diskspace

        Parameters
        ----------
        origin : Origin()
            the origin to be saved
        idx : int or None
            if idx is not None the index will be used for saving in the storage. This might overwrite already existing trajectories!

        Notes
        -----
        This also saves all contained frames in the origin if not done yet.
        A single Origin object can only be saved once!
        """

        idx = super(OriginStorage, self).index(origin, idx)

        if idx is not None:
            storage = self.storage
            storage.variables['origin_name'][idx] = origin.name
            values = []
            for frame_idx, trajectory in enumerate(origin.inputs):
                storage.trajectory.save(trajectory)
                values.append(trajectory.idx[storage])

            self.save_mixed('origin_input', idx, values)

            self.storage.variables['origin_name'][idx] = origin.name
            self.storage.variables['origin_success'][idx] = int(origin.success)
            self.storage.variables['origin_acceptance'][idx] = origin.acceptance
            self.storage.variables['origin_mover_name'][idx] = origin.mover.name

            self.storage.trajectory.save(origin.final)
            self.save_object('origin_final', idx, origin.final)

            self.storage.trajectory.save(origin.result)
            self.save_object('origin_result', idx, origin.result)

            self.storage.ensemble.save(origin.ensemble)
            self.save_object('origin_ensemble', idx, origin.ensemble)

    def load(self, idx, momentum = True):
        '''
        Return a origin from the storage

        Parameters
        ----------
        idx : int
            index of the origin (counts from 1)

        Returns
        -------
        origin : Origin
            the origin
        '''
        name = self.storage.variables['origin_name'][idx]
        trajectories_idcs = self.storage.variables['origin_input'][self.slice('origin_input', idx)]
        success = self.storage.variables['origin_success'][idx]
        acceptance = self.storage.variables['origin_acceptance'][idx]
        mover_name = self.storage.variables['origin_mover_name'][idx]
        final_idx = self.storage.variables['origin_final_idx'][idx]
        result_idx = self.storage.variables['origin_result_idx'][idx]
        ensemble_idx = self.storage.variables['origin_ensemble_idx'][idx]

        obj = Origin(
            name=name,
            final=self.storage.trajectory.load(final_idx, lazy=True),
            result=self.storage.trajectory.load(result_idx, lazy=True),
            inputs=[self.storage.trajectory.load(idx, lazy=True) for idx in trajectories_idcs],
            mover=mover_name,
            acceptance=float(acceptance),
            success=bool(success),
            options=None,
            ensemble=self.storage.ensemble.load(ensemble_idx)
        )

        obj.idx[self.storage] = idx

        return obj

    def _init(self):
        """
        Initialize the associated storage to allow for origin storage

        """
        super(OriginStorage, self)._init()

        # New short-hand definition
        self.init_variable('origin_name', 'str')
        self.init_variable('origin_success', 'b')
        self.init_variable('origin_final_idx', 'u4')
        self.init_variable('origin_result_idx', 'u4')
        self.init_variable('origin_acceptance', 'f')
        self.init_variable('origin_mover_name', 'str')
        self.init_mixed_length('origin_input')
        self.init_variable('origin_input', 'u4', 'origin_input')
        self.init_variable('origin_ensemble_idx', 'u4')
