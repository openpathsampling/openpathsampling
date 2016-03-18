
import random
import logging
import abc

import numpy as np

import openpathsampling as paths
from openpathsampling.netcdfplus import StorableNamedObject, StorableObject

logger = logging.getLogger(__name__)

class SnapshotModifier(StorableNamedObject):
    """Abstract class for snapshot modification.

    In general, a snapshot modifer will take a snapshot and return a
    modified version of that snapshot. These are useful for, e.g., two-way
    shooting or for committor analysis.

    Attributes
    ----------
    subset_mask : list of int or None
        the subset to use (default None, meaning no subset). The values
        select along the first axis of the input array. For example, in a
        typical shape=(n_atoms, 3) array, this will pick the atoms.

    Note
    ----
    It is tempting to use the various indexing tricks in numpy to get a view
    on the data, instead of using this subset_mask. However, this can be
    dangerous: note that fancy indexing copies the data, whereas normal
    indexing returns a proper view. See http://stackoverflow.com/a/4371049.
    Rather than play with fire here, we'll just do something
    straightforward. We can try something more clever in the future.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, subset_mask=None):
        self.subset_mask = subset_mask

    def extract_subset(self, full_array):
        return [full_array[i] for i in self.subset_mask]

    def apply_to_subset(self, full_array, modified):
        for (i, val) in zip(self.subset_mask, modified):
            full_array[i] = val
        return full_array

    @abc.abstractmethod
    def __call__(self, snapshot): # pragma: no cover
        pass

class NoModification(SnapshotModifier):
    """Modifier with no change: returns a copy of the snapshot."""
    def __call__(self, snapshot):
        return snapshot.copy()

class RandomVelocities(SnapshotModifier):
    """Randomize velocities according to the Boltzmann distribution."""
    def __init__(self, subset_mask=None, beta=None):
        pass

    def __call__(self, snapshot):
        pass
