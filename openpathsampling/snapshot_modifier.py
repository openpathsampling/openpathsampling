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
        """Extracts elements from full_array according to self.subset_mask

        Parameters
        ----------
        full_array : list-like
            the input array

        Returns
        -------
        list
            the elements of full_array which are selected by
            self.subset_mask, or full_array if self.subset_mask is None
        """
        if self.subset_mask is None:
            return full_array
        else:
            return [full_array[i] for i in self.subset_mask]

    def apply_to_subset(self, full_array, modified):
        """Replaces elements of full_array according to self.subset_mask

        This returns the full_array, but the modification is done in-place.

        Parameters
        ----------
        full_array : list-like
            array to modify
        modified : list-like
            array containing len(self.subset_mask) elements which will
            replace those in `full_array`

        Returns
        -------
        full_array : list-like
            modified version of the input array, where the elements
            specified by self.subset_mask have been replaced    
        """
        subset_mask = self.subset_mask
        if self.subset_mask is None:
            subset_mask = range(len(full_array))
        for (i, val) in zip(subset_mask, modified):
            full_array[i] = val
        return full_array

    @abc.abstractmethod
    def __call__(self, snapshot): 
        raise NotImplementedError

class NoModification(SnapshotModifier):
    """Modifier with no change: returns a copy of the snapshot."""
    def __call__(self, snapshot):
        return snapshot.copy()

class RandomVelocities(SnapshotModifier):
    """Randomize velocities according to the Boltzmann distribution.
    
    Parameters
    ----------
    beta : float
        Inverse temperature in units of inverse velocity squared [also known
        as 1.0/(k_B * T)].
    """
    def __init__(self, subset_mask=None, beta=None):
        super(RandomVelocities, self).__init__(subset_mask)
        self.beta = beta

    def __call__(self, snapshot):
        new_snap = snapshot.copy()

        vel_subset = self.extract_subset(new_snap.velocities)
        masses = self.extract_subset(snapshot.topology.masses)
        n_spatial = len(vel_subset[0])
        n_atoms = len(vel_subset)
        for atom_i in range(n_atoms):
            sigma = np.sqrt(1.0 / (self.beta * masses[atom_i]))
            vel_subset[atom_i] = sigma * np.random.normal(size=n_spatial)
        self.apply_to_subset(new_snap.velocities, vel_subset)
        return new_snap
