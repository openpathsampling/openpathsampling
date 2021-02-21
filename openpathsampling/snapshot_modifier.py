import random
import logging
import copy
import abc
import functools

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
        the subset to use (default None, meaning use all). The values
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

    def extract_subset(self, full_array, subset=None):
        """Extracts elements from full_array according to self.subset_mask

        Note that, if ``subset`` and ``self.subset`` are None, this returns
        ``full_array``. If you intend to modify the object returned by this
        functions, you should ensure that your input ``full_array`` is a
        copy of any orginal immutable data.

        Parameters
        ----------
        full_array : list-like
            the input array
        subset : list of int or None
            the subset to use; see ``SnapshotModifier.subset_mask``. Default
            (None) uses the value of ``self.subset_mask``.

        Returns
        -------
        list-like
            the elements of full_array which are selected by
            self.subset_mask, or full_array if subset and self.subset_mask
            are None
        """
        if subset is None:
            subset = self.subset_mask

        if subset is None:
            return full_array
        else:
            return [full_array[i] for i in subset]

    def apply_to_subset(self, full_array, modified, subset_mask=None):
        """Replaces elements of full_array according to the subset mask

        This returns the full_array, but the modification is done in-place.

        Parameters
        ----------
        full_array : list-like
            array to modify
        modified : list-like
            array containing len(self.subset_mask) elements which will
            replace those in `full_array`
        subset_mask : list of int or None
            the subset to use; see ``SnapshotModifier.subset_mask``. Default
            (None) uses the value of ``self.subset_mask``.

        Returns
        -------
        full_array : list-like
            modified version of the input array, where the elements
            specified by self.subset_mask have been replaced
        """
        if subset_mask is None:
            subset_mask = self.subset_mask

        if subset_mask is None:
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

    Notes
    -----
    This modifier will only work with snapshots that have the ``velocities``
    feature and the ``masses`` feature. Furthermore, the units have to be
    such that the input ``beta`` and the features `masses` and `velocities`
    are all in the same unit system. In particular, ``1.0 / beta * masses``
    must be in units of ``velocity**2``.

    For the OpenMMEngine, for example (after ``from simtk import unit as
    u``), the ``beta`` parameter for 300 K would be created with

    .. code-block:: python

        beta = 1.0 / (300.0 * u.kelvin * u.BOLTZMANN_CONSTANT_kB)

    Parameters
    ----------
    beta : float or simtk.unit.Quantity
        inverse temperature (including kB) for the distribution
    engine : :class:`.DynamicsEngine` or None
        engine to be used for constraints; if None, use the snapshot's
        engine
    subset_mask : list of int or None
        the subset to use (default None, meaning use all). The values
        select along the first axis of the input array. For example, in a
        typical shape=(n_atoms, 3) array, this will pick the atoms.
    """
    def __init__(self, beta=None, engine=None, subset_mask=None):
        super(RandomVelocities, self).__init__(subset_mask)
        self.beta = beta
        self.engine = engine

    def _default_random_velocities(self, snapshot, beta, subset):
        if beta is None:
            raise RuntimeError("Engine can't use RandomVelocities")

        # raises AttributeError is snapshot doesn't support velocities
        velocities = copy.copy(snapshot.velocities)  # copy.copy for units
        vel_subset = self.extract_subset(velocities, subset)

        # raises AttributeError if snapshot doesn't support masses feature
        all_masses = snapshot.masses
        masses = self.extract_subset(all_masses, subset)

        n_spatial = len(vel_subset[0])
        n_atoms = len(vel_subset)
        for atom_i in range(n_atoms):
            radicand = 1.0 / (self.beta * masses[atom_i])
            try:  # preference that masses is quantity or np.array
                sigma = radicand.sqrt()
            except AttributeError:  # if masses regular list
                sigma = np.sqrt(radicand)
            vel_subset[atom_i] = sigma * np.random.normal(size=n_spatial)

        self.apply_to_subset(velocities, vel_subset, subset)
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

    def __call__(self, snapshot):
        # default value; we'll override if needed
        try:
            beta = self.beta if self.beta is not None else self.engine.beta
        except AttributeError:
            beta = None

        make_snapshot = functools.partial(
            self._default_random_velocities,
            beta=beta,
            subset=self.subset_mask
        )
        if self.engine:
            try:
                make_snapshot = functools.partial(
                    self.engine.randomize_velocities,
                    beta=beta,
                    subset=self.subset_mask
                )
            except AttributeError:
                pass  # use default

        new_snap = make_snapshot(snapshot)
        return new_snap


class GeneralizedDirectionModifier(SnapshotModifier):
    """
    Snapshot modifier which changes velocity direction with constant energy.

    Abstract class with core implementation. The implementation here
    requires takes `delta_v` to be the standard deviation of a Gaussian
    distribution of velocity displacements. The velocities are modified
    according to that displacement, then renormalized so that each atom
    still has the same total momentum as before (although the direction can
    be changed).

    Parameters
    ----------
    delta_v : float or array-like of float
        velocity change parameter, as the width of a Gaussian from which
        each degree of freedom's velocity change is chosen. Can be a list
        with length ``n_atoms`` (mapping to each atom); length equal to the
        subset mask (mapping to each element of the subset mask); or a float
        (mapping the same value to all atoms).
    subset_mask : list of int or None
        the subset to use (default None, meaning use all). The values
        select along the first axis of the input array. For example, in a
        typical shape=(n_atoms, 3) array, this will pick the atoms.
    remove_linear_momentum : bool
        whether the total linear momentum should be removed from the system,
        default is True

    See Also
    --------
    VelocityDirectionModifier
    SingleAtomVelocityDirectionModifier
    """
    def __init__(self, delta_v, subset_mask=None,
                 remove_linear_momentum=True, engine=None):
        super(GeneralizedDirectionModifier, self).__init__(subset_mask)
        self.delta_v = delta_v
        self.remove_linear_momentum = remove_linear_momentum
        self.engine = engine

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
            has_constraints = self.engine.has_constraints()
        except AttributeError:
            pass
        else:
            if has_constraints:
                raise RuntimeError("Cannot use this snapshot modifier with "
                                   "an engine that has constraints.")

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
        """
        Generate the list of velocity delta widths.

        Parameters
        ----------
        n_atoms : int
            number of total atoms
        n_subset_atoms : int
            number of atoms in the subset to be (possibly) changed

        Returns
        -------
        list, length n_subset_atoms
            list of velocity deltas associated with each atom from the
            subset
        """
        try:
            dv_widths = list(self.delta_v)
        except TypeError:
            dv_widths = [self.delta_v] * n_subset_atoms

        if len(dv_widths) == n_atoms:
            dv_widths = self.extract_subset(dv_widths)

        # assert len(dv_widths) == n_subset_atoms
        return dv_widths

    @staticmethod
    def _remove_linear_momentum(velocities, masses):
        """
        Remove COM motion.

        Parameters
        ----------
        velocities : array-like, shape (n_atoms, n_spatial)
            input velocities (after snapshot change)
        masses : array-like, shape (n_atoms,)
            masses of each atom

        Returns
        -------
        array-like
            velocities adjusted to have 0 linear momentum
        """
        # TODO: initially, maybe see if there's an internal motion remover
        # to do most of this? and get KE from a snapshot feature?
        n_atoms = len(masses)
        inv_masses = 1.0 / masses
        momenta = velocities * masses[:, np.newaxis]
        total_momenta = sum(momenta, 0*momenta[0])
        remove_momenta = total_momenta / n_atoms
        remove_velocities = inv_masses[:, np.newaxis] * remove_momenta

        velocities -= remove_velocities
        return velocities

    @staticmethod
    def _rescale_kinetic_energy(velocities, masses, double_KE):
        """
        Rescale KE to desired value.

        Parameters
        ----------
        velocities : array-like, shape (n_atoms, n_spatial)
            input velocities (after snapshot change)
        masses : array-like, shape (n_atoms,)
            masses of each atom
        double_KE : float or unitted Quantity
            the desired kinetic energy multiplied by 2.0 (because to avoid
            needing to multiple by 1/2 internally)

        Returns
        -------
        array-like
            velocities adjusted to have the desired kinetic energy
        """
        # from here, we're doing the KE rescaling
        # can't just use the dot product because of simtk.units
        momenta = velocities * masses[:, np.newaxis]
        dof_ke = momenta * velocities
        zero_energy = 0 * dof_ke[0][0]
        new_ke = sum(sum(dof_ke, zero_energy), zero_energy)

        rescale_factor = np.sqrt(double_KE / new_ke)
        velocities *= rescale_factor

        return velocities

    def __call__(self, snapshot):
        """
        Primary call function to be used by subclasses.

        Uses the :meth:`._select_atoms_to_modify` method to determine
        exactly which atoms from the subset should be modified. This method
        must be written in subclasses.

        Parameters
        ----------
        snapshot : :class:`.Snapshot`
            initial snapshot

        Returns
        -------
        :class:`.Snapshot`
            modified snapshot
        """
        self._verify_snapshot(snapshot)
        velocities = copy.copy(snapshot.velocities)
        vel_subset = self.extract_subset(velocities)

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

        # calculate the total KE so we can preserve it
        masses = snapshot.masses
        momenta = velocities * masses[:, np.newaxis]
        dof_double_KE = momenta * velocities
        zero_energy = 0 * dof_double_KE[0][0]
        double_KE = sum(sum(dof_double_KE, zero_energy), zero_energy)

        if self.remove_linear_momentum:
            velocities = self._remove_linear_momentum(velocities, masses)

        self._rescale_kinetic_energy(velocities, masses, double_KE)

        new_snap = snapshot.copy_with_replacement(velocities=velocities)

        # NOTE: no constraint correction here! constraints are not allowed!
        return new_snap

class VelocityDirectionModifier(GeneralizedDirectionModifier):
    """
    Randomly change all the velocities (in the subset masked, at least)
    according to the given delta_v.

    Parameters
    ----------
    delta_v : float or array-like of float
        velocity change parameter, as the width of a Gaussian from which
        each degree of freedom's velocity change is chosen. Can be a list
        with length ``n_atoms`` (mapping to each atom); length equal to the
        subset mask (mapping to each element of the subset mask); or a float
        (mapping the same value to all atoms).
    subset_mask : list of int or None
        the subset to use (default None, meaning use all). The values
        select along the first axis of the input array. For example, in a
        typical shape=(n_atoms, 3) array, this will pick the atoms.
    rescale_linear_momenta : bool
        whether the total linear momentum should be removed from the system,
        default is True

    See Also
    --------
    GeneralizedDirectionModifier
    SingleAtomVelocityDirectionModifier
    """
    def _select_atoms_to_modify(self, n_subset_atoms):
        return list(range(n_subset_atoms))

class SingleAtomVelocityDirectionModifier(GeneralizedDirectionModifier):
    """
    Change a single atom according to the ``delta_v``.

    Note that even though this only uses the ``delta_v`` to change one atom,
    the velocities of other atoms will also be changed in order to preserve
    the initial kinetic energy and possibly also to remove linear momenta.

    Parameters
    ----------
    delta_v : float or array-like of float
        velocity change parameter, as the width of a Gaussian from which
        each degree of freedom's velocity change is chosen. Can be a list
        with length ``n_atoms`` (mapping to each atom); length equal to the
        subset mask (mapping to each element of the subset mask); or a float
        (mapping the same value to all atoms).
    subset_mask : list of int or None
        the subset to use (default None, meaning use all). The values
        select along the first axis of the input array. For example, in a
        typical shape=(n_atoms, 3) array, this will pick the atoms.
    rescale_linear_momenta : bool
        whether the total linear momentum should be removed from the system,
        default is True

    See Also
    --------
    GeneralizedDirectionModifier
    VelocityDirectionModifier
    """
    def _select_atoms_to_modify(self, n_subset_atoms):
        return [np.random.choice(range(n_subset_atoms))]

