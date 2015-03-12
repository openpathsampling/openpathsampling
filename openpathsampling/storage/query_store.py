class QueryStore():

    def __init__(self, storage):
        self.storage = storage


    def orderparameter_trajectory(self, orderparameter, ensemble=None, replica=None, step=None, trial=False):

        storage = self.storage

        ens_idx = ensemble.idx[storage]

        output = []

        op_dict = orderparameter.storage_caches[storage]

        for sset_id in range(len(storage.sampleset)):
            if step is not None and sset_id != step:
                continue

            sample_idxs = storage.variables['sampleset_sample_idx'][sset_id].tolist()

            ensemble_idxs = storage.variables['sample_ensemble_idx'][sample_idxs].tolist()
            replica_idxs = storage.variables['sample_replica_idx'][sample_idxs].tolist()
            traj_idx = storage.variables['sample_trajectory_idx'][sample_idxs].tolist()

            for no, sample_idx in enumerate(sample_idxs):
                if ensemble is not None and ens_idx != ensemble_idxs[no]:
                    continue
                if replica is not None and replica != replica_idxs[no]:
                    continue

                snap_idxs = storage.variables['trajectory_snapshot_idx'][traj_idx[no]]

                output.append([ op_dict[idx] for idx in snap_idxs ])

        return output