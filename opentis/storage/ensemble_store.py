from object_storage import ObjectStorage
from opentis.ensemble import Ensemble, LoadedEnsemble
from wrapper import loadcache, savecache, savenamed

class EnsembleStorage(ObjectStorage):

    def __init__(self, storage):
        super(EnsembleStorage, self).__init__(storage, Ensemble, named=True)

    @loadcache
    def load(self, idx):
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

        name = self.load_variable('ensemble_name', idx)
        description = self.load_variable('ensemble_str', idx)

        obj = LoadedEnsemble(name=name, description=description)
        return obj


    @savenamed
    @savecache
    def save(self, ensemble, idx):
        '''
        Returns an object from the storage. Needs to be implemented from the specific storage class.
        '''

        if self.named and hasattr(ensemble, 'name'):
            self.save_variable(self.db + '_name', idx, ensemble.name)

        self.save_variable('ensemble_str', idx, str(ensemble))


    def _init(self):
        """
        Initialize the associated storage to allow for ensemble storage

        """
        super(EnsembleStorage, self)._init()

        self.init_variable('ensemble_str', 'str', description="ensemble_str[ensemble] is the description string of ensemble 'ensemble'.")
