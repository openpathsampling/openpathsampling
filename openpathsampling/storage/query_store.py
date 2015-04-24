class QueryStore():

    def __init__(self, storage):
        self.storage = storage


    def trajectory_collectivevariable(self, collectivevariable, ensemble=None, replica=None, step=None):
        """
        Return list of collectivevariables fast for specific sets of samples
        samples can be all samples found in specific or all sample_set and filter
        these by ensemble and/or replica.
        Parameters
        ----------
        collectivevariable : paths.CollectiveVariable()
            the collectivevariable from which the values should be extracted
        ensemble : paths.Ensemble or None
            if not None only samples from the specific ensemble are used.
            For `None` (default) all ensembles are considered
        replica : int or None
            if not None only samples from the specific replica ID are used.
            For `None` (default) all replica IDs are considered.
        step : int or None
            if not None only samples from the specific step are used.
            For `None` (default) all sample_set steps are considered.
        Returns
        -------
        list of list of float
            Returns for each sample a list of floats which represent the
            collectivevariable values of the trajectory of the samples
        """

        storage = self.storage

        if ensemble is not None:
            ens_idx = ensemble.idx[storage]

        output = []

        op_dict = storage.cvs.get_list_value(collectivevariable, slice(None,None))

        for sset_id in range(len(storage.sample_set)):
            if step is not None and sset_id != step:
                continue

            sample_idxs = storage.variables['sampleset_sample_idx'][sset_id].tolist()
            ensemble_idxs = storage.variables['sample_ensemble_idx'][sample_idxs].tolist()
            replica_idxs = storage.variables['sample_replica'][sample_idxs].tolist()
            traj_idx = storage.variables['sample_trajectory_idx'][sample_idxs].tolist()

            for no, sample_idx in enumerate(sample_idxs):
                if ensemble is not None and ens_idx != ensemble_idxs[no]:
                    continue
                if replica is not None and replica != replica_idxs[no]:
                    continue

                snap_idxs = storage.variables['trajectory_snapshot_idx'][traj_idx[no]]

                output.append([ op_dict[idx] for idx in snap_idxs ])

        return output

    def trajectory_length(self, ensemble=None, replica=None, step=None):
        """
        Return list of trajectory lengths fast for specific sets of samples
        samples can be all samples found in specific or all sample_set and filter
        these by ensemble and/or replica.
        Parameters
        ----------
        ensemble : paths.Ensemble or None
            if not None only samples from the specific ensemble are used.
            For `None` (default) all ensembles are considered
        replica : int or None
            if not None only samples from the specific replica ID are used.
            For `None` (default) all replica IDs are considered.
        step : int or None
            if not None only samples from the specific step are used.
            For `None` (default) all sample_set steps are considered.
        Returns
        -------
        list of int
            Returns for each sample the number of snapshots in it
        """
        storage = self.storage

        if ensemble is not None:
            ens_idx = ensemble.idx[storage]

        output = []

        for sset_id in range(len(storage.sample_set)):
            if step is not None and sset_id != step:
                continue

            sample_idxs = storage.variables['sampleset_sample_idx'][sset_id].tolist()
            ensemble_idxs = storage.variables['sample_ensemble_idx'][sample_idxs].tolist()
            replica_idxs = storage.variables['sample_replica'][sample_idxs].tolist()
            traj_idx = storage.variables['sample_trajectory_idx'][sample_idxs].tolist()

            for no, sample_idx in enumerate(sample_idxs):
                if ensemble is not None and ens_idx != ensemble_idxs[no]:
                    continue
                if replica is not None and replica != replica_idxs[no]:
                    continue

                snap_idxs = storage.variables['trajectory_snapshot_idx'][traj_idx[no]]

                output.append(len(snap_idxs))

        return output