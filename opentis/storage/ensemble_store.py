from object_storage import ObjectStorage
from opentis.ensemble import Ensemble, LoadedEnsemble
from wrapper import loadcache

class EnsembleStorage(ObjectStorage):

    def __init__(self, storage):
        super(EnsembleStorage, self).__init__(storage, Ensemble, named=True)

    def save(self, ensemble, idx=None):
        """
        Add the current state of the ensemble in the database. If nothing has changed then the ensemble gets stored using the same snapshots as before. Saving lots of diskspace

        Parameters
        ----------
        ensemble : Ensemble()
            the ensemble to be saved
        idx : int or None
            if idx is not None the index will be used for saving in the storage. This might overwrite already existing trajectories!

        Notes
        -----
        This also saves all contained frames in the ensemble if not done yet.
        A single Ensemble object can only be saved once!
        """

        if idx is None:
            find_idx = self.find_by_name(str(ensemble))
            if find_idx is not None:
                # found and does not need to be saved, but we will let this ensemble point to the storage
                # in case we want to save and need the idx
                ensemble.idx[self.storage] = find_idx
                self.cache[find_idx] = ensemble
                # TODO: We might check if the string representation agrees and throw an exception otherwise
                return

        idx = super(EnsembleStorage, self).index(ensemble, idx)

        if idx is not None:
            storage = self.storage
            storage.variables['ensemble_name'][idx] = ensemble.name
            storage.variables['ensemble_str'][idx] = str(ensemble)

            ensemble.idx[storage] = idx

        return

    def find_by_name(self, needle):
        all_names = self.storage.variables['ensemble_str'][:]
        for idx, name in enumerate(all_names):
            if name == needle:
                return idx

        return None

    @loadcache
    def load(self, idx, momentum = True):
        '''
        Return a ensemble from the storage

        Parameters
        ----------
        idx : int
            index of the ensemble (counts from 1)

        Returns
        -------
        ensemble : Ensemble
            the ensemble
        '''

        name = self.storage.variables['ensemble_name'][int(idx)]
        description = self.storage.variables['ensemble_str'][int(idx)]

        obj = LoadedEnsemble(name=name, description=description)
        obj.idx[self.storage] = idx

        return obj


    def _init(self):
        """
        Initialize the associated storage to allow for ensemble storage

        """
        super(EnsembleStorage, self)._init()

        self.init_variable('ensemble_str', 'str', description="ensemble_str[ensemble] is the description string of ensemble 'ensemble'.")
