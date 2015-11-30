from openpathsampling.sample import SampleSet, Sample
from openpathsampling.netcdfplus import VariableStore


class SampleStore(VariableStore):
    def __init__(self):
        super(SampleStore, self).__init__(
            Sample,
            ['trajectory', 'ensemble', 'replica', 'parent', 'details', 'bias', 'mover']
        )

    def by_ensemble(self, ensemble):
        return [sample for sample in self.iterator() if sample.ensemble == ensemble]

    def _init(self):
        super(SampleStore, self)._init()

        # New short-hand definition
        self.create_variable('trajectory', 'obj.trajectories')
        self.create_variable('ensemble', 'obj.ensembles')
        self.create_variable('replica', 'int')
        self.create_variable('parent', 'lazyobj.samples')
        self.create_variable('details', 'lazyobj.details')
        self.create_variable('bias', 'float')
        self.create_variable('mover', 'obj.pathmovers')


class SampleSetStore(VariableStore):
    def __init__(self):
        super(SampleSetStore, self).__init__(
            SampleSet,
            ['samples', 'movepath']
        )


    def sample_indices(self, idx):
        """
        Load sample indices for sample_set with ID 'idx' from the storage

        Parameters
        ----------
        idx : int
            ID of the sample_set

        Returns
        -------
        list of int
            list of sample indices
        """

        return self.variables['samples'][idx].tolist()

    def _init(self, units=None):
        """
        Initialize the associated storage to allow for sampleset storage

        """
        super(SampleSetStore, self)._init()

        self.create_variable('samples', 'obj.samples',
                           dimensions='...',
                           description="sampleset[sampleset][frame] is the sample index " +
                                       "(0..nspanshots-1) of frame 'frame' of sampleset 'sampleset'.",
                           chunksizes=(1024,)
                           )

        self.create_variable('movepath', 'lazyobj.pathmovechanges')
