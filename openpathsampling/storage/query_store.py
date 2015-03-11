import numpy as np

from object_storage import ObjectStore
from openpathsampling.trajectory import Trajectory


class QueryStore():

    def __init__(self, storage):
        self.storage = storage

    def orderparameter_samples(self, ensemble, op):
        storage = self.storage
        ensemble_idx = ensemble.idx[storage]
        op_idx = op.idx[storage]