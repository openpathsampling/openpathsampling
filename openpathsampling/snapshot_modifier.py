import random
import logging
import copy
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
        super(SnapshotModifier, self).__init__()
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

    Note
    ----
    This modifier will only work with snapshots that have the `velocities`
    feature and the `masses` feature. Furthermore, the units have to be such
    that the input `beta` and the features `masses` and `velocities` are all
    in the same unit system. In particular, `1.0 / beta * masses` must be in
    units of `velocity**2`.

    Parameters
    ----------
    beta : float
        inverse temperature (in units of kB) for the distribution
    engine : :class:`.DynamicsEngine` or None
        engine to be used for constraints; if None, use the snapshot's
        engine
    subset_mask : list of int or None
        the subset to use (default None, meaning no subset). The values
        select along the first axis of the input array. For example, in a
        typical shape=(n_atoms, 3) array, this will pick the atoms.
    """
    def __init__(self, beta, engine=None, subset_mask=None):
        super(RandomVelocities, self).__init__(subset_mask)
        self.beta = beta
        self.engine = engine

    def __call__(self, snapshot):
        # raises AttributeError is snapshot doesn't support velocities
        velocities = copy.copy(snapshot.velocities)  # copy.copy for units
        vel_subset = self.extract_subset(velocities)

        # raises AttributeError if snapshot doesn't support masses feature
        all_masses = snapshot.masses
        masses = self.extract_subset(all_masses)

        n_spatial = len(vel_subset[0])
        n_atoms = len(vel_subset)
        for atom_i in range(n_atoms):
            sigma = np.sqrt(1.0 / (self.beta * masses[atom_i]))
            vel_subset[atom_i] = sigma * np.random.normal(size=n_spatial)

        self.apply_to_subset(velocities, vel_subset)
        new_snap = snapshot.copy_with_replacement(velocities=velocities)

        # applying constraints, if they exist
        if self.engine is None:
            engine = new_snap.engine
        else:
            engine = self.engine

        try:
            apply_constraints = engine.apply_constraints
        except AttributeError:
            pass  # fine if there isn't one
        else:
            new_snap = apply_constraints(new_snap)

        return new_snap

class GeneralizedDirectionModifier(SnapshotModifier):
    """
    Snapshot modifier which changes momentum direction with constant energy.
    """
    def __init__(self, delta_v, subset_mask=None):
        super(GeneralizedDirectionModifier, self).__init__(subset_mask)
        self.delta_v = delta_v

    def _verify_snapshot(self, snapshot):
        """
        Verifies that a snapshot has the right number of degrees of freedom.

        The approach implemented in this will not satisfy detailed balance
        if there are constraints on the atoms that are changing. This is a
        way of checking for that.

        Parameters
        ----------
        snapshot : :class:`.Snapshot`
            input snapshot to check for validity
        """
        try:
            box_vectors = snapshot.box_vectors
        except AttributeError:
            box_vectors = None

        try:
            n_dofs = snapshot.n_degrees_of_freedom
        except AttributeError:
            raise RuntimeError("Snapshot missing n_degrees_of_freedom. "
                               + "Can't use this snapshot modifier.")

        # all engines should have n_spatial and n_atoms
        n_spatial = snapshot.engine.n_spatial
        n_atoms = snapshot.engine.n_atoms

        # NOTE: none of our engines currently explicitly remove angular
        # (but isn't angular momentum impossible in periodic condensed
        # phase)
        remove_angular = 0  # if box_vectors is None else 1
        remove_linear = 1 if n_atoms != 1 else 0
        if remove_linear:
            try:
                ignore = snapshot.engine.ignore_linear_momentum
            except AttributeError:
                ignore = False
            remove_linear = 0 if ignore else remove_linear

        n_motion_removers = n_spatial * (remove_linear + remove_angular)

        n_dofs_required = n_spatial * n_atoms - n_motion_removers

        if n_dofs != n_dofs_required:
            raise RuntimeError("Snapshot has " + str(n_dofs)
                               + " degrees of freedom, not " 
                               + str(n_dofs_required) + ". "
                               + "Are there constraints? Constraints can't"
                               + " be used with this modifier.")

    def _select_atoms_to_modify(self, n_subset_atoms):
        raise NotImplementedError

    def _dv_widths(self, n_atoms, n_subset_atoms):
        try:
            dv_widths = list(self.delta_v)
        except TypeError:
            dv_widths = [self.delta_v] * n_subset_atoms

        if len(dv_widths) == n_atoms:
            dv_widths = self.extract_subset(dv_widths)

        # assert len(dv_widths) == n_subset_atoms
        return dv_widths

    def __call__(self, snapshot):
        self._verify_snapshot(snapshot)
        velocities = copy.copy(snapshot.velocities)
        vel_subset = self.extract_subset(velocities)
        masses = self.extract_subset(snapshot.masses)

        atoms_to_change = self._select_atoms_to_modify(len(vel_subset))
        dv_widths = self._dv_widths(n_atoms=len(velocities),
                                    n_subset_atoms=len(vel_subset))

        # zero_with_units is a hack for using units (or not) in `sum`
        zero_with_units = (vel_subset[0][0] - vel_subset[0][0])**2
        for atom_i in atoms_to_change:
            initial_sum_sq_vel = sum([v**2 for v in vel_subset[atom_i]],
                                     zero_with_units)
            randoms = np.random.normal(size=len(vel_subset[atom_i]))
            delta_v = dv_widths[atom_i] * randoms
            vel_subset[atom_i] += delta_v
            final_sum_sq_vel = sum([v**2 for v in vel_subset[atom_i]],
                                   zero_with_units)
            rescale_factor = np.sqrt(initial_sum_sq_vel / final_sum_sq_vel)
            vel_subset[atom_i] *= rescale_factor

        self.apply_to_subset(velocities, vel_subset)
        new_snap = snapshot.copy_with_replacement(velocities=velocities)

        # NOTE: no constraint correction here! constraints are not allowed!
        return new_snap

class VelocityDirectionModifier(GeneralizedDirectionModifier):
    def _select_atoms_to_modify(self, n_subset_atoms):
        return range(n_subset_atoms)

class SingleAtomVelocityDirectionModifier(GeneralizedDirectionModifier):
    def _select_atoms_to_modify(self, n_subset_atoms):
        return [np.random.choice(range(n_subset_atoms))]

