from object_storage import ObjectStorage

class Sample(object):
    def __init__(self, trajectory, ensemble, origin):
        self.trajectory = trajectory
        self.ensemble = ensemble
        self.origin = origin
        pass


class SampleStorage(ObjectStorage):

    def __init__(self, storage):
        super(SampleStorage, self).__init__(storage, Sample)

    def save(self, sample, idx=None):
        """
        Add the current state of the sample in the database. If nothing has changed then the sample gets stored using the same snapshots as before. Saving lots of diskspace

        Parameters
        ----------
        sample : Sample()
            the sample to be saved
        idx : int or None
            if idx is not None the index will be used for saving in the storage. This might overwrite already existing trajectories!

        Notes
        -----
        This also saves all contained frames in the sample if not done yet.
        A single Sample object can only be saved once!
        """

        idx = super(SampleStorage, self).index(sample, idx)

        if idx is not None:
            storage = self.storage

            storage.trajectory.save(sample.trajectory)
            storage.ensemble.save(sample.ensemble)
            storage.origin.save(sample.origin)

            storage.variables['sample_trajectory_idx'][idx] = sample.trajectory.idx[self.storage]
            storage.variables['sample_ensemble_idx'][idx] = sample.ensemble.idx[self.storage]
            storage.variables['sample_origin_idx'][idx] = sample.origin.idx[self.storage]

        return

    def load(self, idx, momentum = True):
        '''
        Return a sample from the storage

        Parameters
        ----------
        idx : int
            index of the sample (counts from 1)

        Returns
        -------
        sample : Sample
            the sample
        '''
        trajectory_idx = int(self.storage.variables['sample_trajectory_idx'][idx])
        ensemble_idx = int(self.storage.variables['sample_ensemble_idx'][idx])
        origin_idx = int(self.storage.variables['sample_origin_idx'][idx])

        obj = Sample(
            trajectory=self.storage.trajectory.load(trajectory_idx),
            ensemble = self.storage.ensemble.load(ensemble_idx),
            origin = self.storage.origin.load(origin_idx)
        )
        obj.idx[self.storage] = idx

        return obj


    def _init(self):
        """
        Initialize the associated storage to allow for sample storage

        """
        super(SampleStorage, self)._init()

        # index associated storage in class variable for all Sample instances to access
        ncfile = self.storage

        # Create variables for trajectories
        ncvar_sample_trajectory_idx     = ncfile.createVariable('sample_trajectory_idx', 'u4', self.idx_dimension)
        ncvar_sample_ensemble_idx       = ncfile.createVariable('sample_ensemble_idx', 'u4', self.idx_dimension)
        ncvar_sample_origin_idx         = ncfile.createVariable('sample_origin_idx', 'u4', self.idx_dimension)


        # Define units for snapshot variables.
        setattr(ncvar_sample_trajectory_idx,      'units', 'none')
        setattr(ncvar_sample_ensemble_idx,        'units', 'none')
        setattr(ncvar_sample_origin_idx,          'units', 'none')

        # Define long (human-readable) names for variables.
        setattr(ncvar_sample_trajectory_idx,    "long_name", "sample[sample] is the trajectory index (0..trajectory-1) of sample 'sample'.")
        setattr(ncvar_sample_ensemble_idx,      "long_name", "sample[sample] is the ensemble index (0..ensemble-1) of sample 'sample'.")
        setattr(ncvar_sample_origin_idx,        "long_name", "sample[sample] is the origin index (0..origin-1) of sample 'sample'.")
