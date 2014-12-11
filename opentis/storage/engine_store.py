from object_storage import ObjectStorage
from wrapper import loadcache, savecache, savenamed
from opentis.dynamics_engine import DynamicsEngine

class DynamicsEngineStorage(ObjectStorage):

    def __init__(self, storage):
        super(DynamicsEngineStorage, self).__init__(storage, DynamicsEngine, named=True)

    @loadcache
    def load(self, idx, lazy=True):
        '''
        Returns an object from the storage. Needs to be implented from the specific storage class.
        '''

        # Create object first to break any unwanted recursion in loading
        engine_type = self.load_variable('dynamicsengine_type', idx)

        # try to create an object of the same type as the original
        engine_class_dict = DynamicsEngine.__descendents__()

        options = self.simplifier.from_json(self.load_variable('dynamicsengine_options', idx))
        engine = engine_class_dict[engine_type](options=options)

        return engine

    @savenamed
    @savecache
    def save(self, engine, idx=None):
        '''
        Returns an object from the storage. Needs to be implented from the specific storage class.
        '''

        self.save_variable('dynamicsengine_type', idx, engine.__class__.__name__)

        if self.named and hasattr(engine, 'name'):
            self.save_variable(self.db + '_name', idx, engine.name)


        self.save_variable('dynamicsengine_options', idx, self.simplifier.to_json(engine.options))

        engine.idx[self.storage] = idx

    def _init(self):
        """
        Initialize the associated storage to allow for ensemble storage

        """
        super(DynamicsEngineStorage, self)._init()

        self.init_variable('dynamicsengine_options', 'str', description="dynamicsengine_options[ensemble] is the json string of the options dict of engine 'ensemble'.")
        self.init_variable('dynamicsengine_type', 'str', description="dynamicsengine_type[ensemble] is the actual class name of engine 'ensemble'.")