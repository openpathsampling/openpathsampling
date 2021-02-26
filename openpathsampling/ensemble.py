"""
Created on 03.09.2014

@author: Jan-Hendrik Prinz, David W.H. Swenson
"""

import abc
import logging
import itertools

from openpathsampling.netcdfplus import StorableNamedObject
import openpathsampling as paths

from future.utils import with_metaclass


logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


# TODO: Make Full and Empty be Singletons to avoid storing them several times!


def join_ensembles(ensemble_list):
    """Join several ensembles using a set theory union.

    Parameters
    ----------
    ensemble_list : list of :class:`.Ensemble`
        list of ensembles to join

    Returns
    -------
    :class:`.Ensemble`
        union of all given ensembles
    """
    ensemble = None
    for ens in ensemble_list:
        if ensemble is None:
            ensemble = ens
        else:
            ensemble = ensemble | ens
    return ensemble


def _get_list_traj(trajectory):
    """Return a list of (proxy) snapshots from either a list or Trajectory

    Parameters
    ----------
    trajectory : :class:`.Trajectory` or :class:`.list`
       trajectory or list to convert into a list of snapshots

    Returns
    -------
    :class:`.list`
       list of (proxy) snapshots

    Note
    ----
    Due to possible UUID space restrictions, we don't want to slice
    Trajectories if we can help it (as this will generate a new Trajectory
    object with its own UUID). This is a convenience function that will either
    turn a Trajectory into a list of (proxy) snapshots or just list to list
    mapping, which you can slice without generating a new UUID

    For a more in-depth discussion, please see:
    https://github.com/openpathsampling/openpathsampling/pull/978
    """
    itraj = getattr(trajectory, 'iter_proxies', trajectory.__iter__)
    return list(itraj())


# note: the cache is not storable, because that would just be silly!
class EnsembleCache(object):
    """Object used by ensembles to enable fast algorithms for basic functions.

    The contents stored in the `can_append`, `can_prepend`, `call`, and
    `check_reverse` dictionaries will depend on the ensemble. Only two of
    these dictionaries should be non-`None` at any time: either the pair
    `call` and `can_append`, or the pair `check_reverse` and `can_prepend`.

    This object also contains basic functions to manage the cache.

    Attributes
    ----------
        start_frame : :class:`openpathsampling.snapshot.Snapshot`
        prev_last_frame : :class:`openpathsampling.snapshot.Snapshot`
        direction : +1 or -1
        contents : dictionary
    """

    def __init__(self, direction=None):
        self.start_frame = None
        self.prev_last_frame = None
        self.prev_last_index = None
        self.last_length = None
        self.direction = direction
        self.contents = {}
        self.trusted = False
        self.debug_enabled = False

    def bad_direction_error(self):
        raise RuntimeError("EnsembleCache.direction = " +
                           str(self.direction) + " invalid.")  # nocover

    # def clear(self):
    #     self.start_frame = None
    #     self.prev_last_frame = None
    #     self.last_length = None
    #     self.contents = {}

    def check(self, trajectory=None, reset=None):
        """Checks and resets (if necessary) the ensemble cache.

        The trajectory is considered trustworthy based on checking several
        factors, compared to the last time the cache was checked. For
        forward caches (direction > 0), these are

        * the first frame has not changed
        * the length is the same, or has changed by 1
        * if length unchanged, the final frame is the same; if length
          changed by 1, the penultimate frame is the old final frame

        Similar rules apply for backward caches (direction < 0), with
        obvious changes of "final" and "first" frames.

        If the trajectory is not trustworthy, we return True (should be
        reset).

        Parameters
        ----------
        trajectory : :class:`.Trajectory`
            the trajectory to test
        reset : bool or None
            force a value for reset. If None, the value is determined based
            on the test criteria.

        Returns
        -------
        bool :
            the value of reset
        """
        if self.debug_enabled:
            logger.debug("Checking cache....")
            # logger.debug("traj " + str([id(s) for s in trajectory]))
            logger.debug("start_frame " + str(id(self.start_frame)))
            logger.debug("prev_last " + str(id(self.prev_last_frame)))
            logger.debug("prev_last_idx " + str(self.prev_last_index))

        if trajectory is not None:
            # this might get a list instead of Trajectory from internal
            # functions
            get_frame = getattr(trajectory, "get_as_proxy",
                                trajectory.__getitem__)

            # if the first frame has changed, we should reset
            if reset is None:
                lentraj = len(trajectory)
                if self.direction > 0:
                    if get_frame(0) != self.start_frame:
                        reset = True
                    else:
                        if lentraj == 1:
                            # makes no difference here; always reset
                            reset = True
                        elif lentraj == self.last_length:
                            reset = (get_frame(-1) != self.prev_last_frame)
                        elif lentraj == self.last_length + 1:
                            reset = (get_frame(-2) != self.prev_last_frame)
                        else:
                            reset = True
                elif self.direction < 0:
                    if get_frame(-1) != self.start_frame:
                        reset = True
                    else:
                        if lentraj == 1:
                            reset = True
                        elif lentraj == self.last_length:
                            reset = (get_frame(0) != self.prev_last_frame)
                        elif lentraj == self.last_length + 1:
                            reset = (get_frame(1) != self.prev_last_frame)
                        else:
                            reset = True
                else:
                    self.bad_direction_error()
        else:
            reset = True

        self.trusted = not reset
        self.last_length = len(trajectory)
        if reset:
            self.debug_enabled = logger.isEnabledFor(logging.DEBUG)
            if self.debug_enabled:
                logger.debug("Resetting cache " + str(self))
            if self.direction > 0:
                # TODO: this can be hit with trajectory is None?
                self.start_frame = get_frame(0)
                self.prev_last_frame = get_frame(-1)
                self.last_length = len(trajectory)
                self.contents = {}
            elif self.direction < 0:
                # TODO: this can be hit with trajectory is None?
                self.start_frame = get_frame(-1)
                self.prev_last_frame = get_frame(0)
                self.last_length = len(trajectory)
                self.contents = {}
            else:
                self.bad_direction_error()
        else:
            self.trusted = True
        # by returning reset, we allow the functions that call this to reset
        # other things as well
        if self.direction > 0:
            # TODO: this can be hit with trajectory is None?
            self.prev_last_frame = get_frame(-1)
            self.prev_last_index = len(trajectory) - 1
        elif self.direction < 0:
            # TODO: this can be hit with trajectory is None?
            self.prev_last_frame = get_frame(0)
            self.prev_last_index = 0
        else:
            self.bad_direction_error()

        return reset


class Ensemble(with_metaclass(abc.ABCMeta, StorableNamedObject)):
    """
    Path ensemble object.

    An Ensemble represents a path ensemble, effectively a set of trajectories.
    Typical set operations are allowed, here: and, or, xor, -(without), ~
    (inverse = all - x)

    Notes
    -----
    Maybe replace - by / to get better notation. So far it has not been used
    """

    #__metaclass__ = abc.ABCMeta

    def __init__(self):
        """
        A path volume defines a set of paths.
        """
        super(Ensemble, self).__init__()
        self._saved_str = None  # cached first time it is requested

    # https://docs.python.org/3/reference/datamodel.html#object.__hash__
    __hash__ = StorableNamedObject.__hash__

    def __eq__(self, other):
        if self is other:
            return True
        return str(self) == str(other)

    def __ne__(self, other):
        return not self == other

    @abc.abstractmethod
    def __call__(self, trajectory, trusted=None, candidate=False):
        """
        Return `True` if the trajectory is part of the path ensemble.

        Parameters
        ----------
        trajectory: :class:`.Trajectory`
            The trajectory to be checked
        trusted : boolean
            For many ensembles, a faster algorithm can be used if we know
            some information about the trajectory with one fewer frames.
            The `trusted` flag tells the ensemble to use such an algorithm.
            This is usually used in combination with an
            :class:`.EnsembleCache` which makes short-cut calculations
            possible.
        """
        return False

    def check_reverse(self, trajectory, trusted=False):
        """
        See __call__; same thing, but potentially in reverse frame order
        """
        return self(trajectory, trusted=False)

    def check(self, trajectory):
        """Alias for __call__"""
        return self(trajectory, trusted=False)

    def trajectory_summary(self, trajectory):
        """
        Return dict with info on how this ensemble "sees" the trajectory.

        Parameters
        ----------
        trajectory : `openpathsampling.Trajectory`
        """
        return {}

    def trajectory_summary_str(self, trajectory):
        """
        Returns a string with the results of the trajectory_summary function.

        Parameters
        ----------
        trajectory : `openpathsampling.Trajectory`
        """
        summ = self.trajectory_summary(trajectory)
        if summ == {}:
            return "No summary available"
        else:
            return str(summ)

    def can_append(self, trajectory, trusted=False):
        """
        Returns true, if the trajectory so far can still be in the ensemble
        if it is appended by a frame. To check, it assumes that the
        trajectory to length L-1 is okay. This is mainly for interactive
        usage, when a trajectory is generated.

        Parameters
        ----------
        trajectory : :class:`openpathsampling.trajectory.Trajectory`
            the actual trajectory to be tested
        trusted : bool
            If trusted=True, some ensembles can be computed more efficiently
            (e.g., by checking only one frame)

        Returns
        -------
        bool
            Returns true or false if using a forward step (extending the
            trajectory forward in time at its end) `trajectory` could  still
            be in the ensemble and thus makes sense to continue a simulation
        """
        return True

    def can_prepend(self, trajectory, trusted=False):
        """
        Returns true, if the trajectory so far can still be in the ensemble
        if it is prepended by a frame. To check, it assumes that the
        trajectory from index 1 is okay. This is mainly for interactive
        usage, when a trajectory is generated using a backward move.

        Parameters
        ----------
        trajectory : :class:`openpathsampling.trajectory.Trajectory`
            the actual trajectory to be tested
        trusted : bool
            If trusted=True, some ensembles can be computed more efficiently
            (e.g., by checking only one frame)

        Returns
        -------
        bool
            Returns true or false if using a backward step (extending the
            trajectory backwards in time at its beginning) `trajectory`
            could  still be in the ensemble and thus makes sense to continue
            a simulation
        """
        return True

    def strict_can_append(self, trajectory, trusted=False):
        """
        Returns true if the trajectory can be the beginning of a trajectory
        in the ensemble.

        Parameters
        ----------
        trajectory : :class:`.Trajectory`
            trajectory to test
        trusted : bool
            If trusted=True, some ensembles can be computed more efficiently
            (e.g., by checking only one frame)

        Returns
        -------
        bool
            True if and only if the given trajectory can be the beginning of
            a trajectory in the ensemble.
        """
        # default behavior is to be the same as can_append
        return self.can_append(trajectory, trusted)

    def strict_can_prepend(self, trajectory, trusted=False):
        """
        Returns true if the trajectory can be the end of a trajectory in the
        ensemble.

        Parameters
        ----------
        trajectory : :class:`.Trajectory`
            trajectory to test
        trusted : bool
            If trusted=True, some ensembles can be computed more efficiently
            (e.g., by checking only one frame)

        Returns
        -------
        bool
            True if and only if the given trajectory can be the end of a
            trajectory in the ensemble.
        """
        # default behavior is to be the same as can_prepend
        return self.can_prepend(trajectory, trusted)

    def iter_valid_slices(
            self,
            trajectory,
            max_length=None,
            min_length=1,
            overlap=1,
            reverse=False
    ):
        """
        Return an iterator over slices of subtrajectories matching the ensemble

        Parameters
        ----------
        trajectory : :class:`openpathsampling.trajectory.Trajectory`
            the actual trajectory to be splitted into ensemble parts
        max_length : int > 0, optional
            if set this determines the maximal size to be tested (is mainly
            used in the recursion)
        min_length : int > 0, optional
            if set this determines the minimal size to be tested (in lazy
            mode might no
        overlap : int >= 0, optional
            determines the allowed overlap of all trajectories to be found.
            A value of x means that two sub-trajectorie can share up to x
            frames at the beginning and x frames at the end.  Default is 1
        reverse : bool
            if `True` this will start searching from the end of the trajectory.
            Otherwise (default) it will start at the beginning.

        Returns
        -------
        list of `slice`
            Returns a list of index-slices for sub-trajectories in
            trajectory that are in the ensemble.
        """
        length = len(trajectory)

        if max_length is None:
            max_length = length

        max_length = min(length, max_length)
        min_length = max(1, min_length)

        logger.debug("Looking for subtrajectories in " + str(trajectory))
        old_tt_len = 0

        if not reverse:
            start = 0
            end = start + min_length

            while start <= length - min_length and end <= length:
                # print start, end
                tt = trajectory[start:end]

                if len(tt) != old_tt_len + 1:
                    can_append_tt = self.strict_can_append(tt)
                else:
                    can_append_tt = self.strict_can_append(tt, trusted=True)
                old_tt_len = len(tt)

                if end < length and can_append_tt:
                    end += 1
                    if end - start > max_length + 1:
                        start += 1
                        end = start + min_length
                else:
                    if end - start <= max_length and self(tt, trusted=False):
                        yield slice(start, end)
                        pad = min(overlap, end - start - 1)
                        start = end - pad
                        if end == length:
                            # This means we have reached the end and should stop
                            # All other possible subtraj can only be contained
                            # in already existing ones
                            start = length
                    elif end - start >= min_length + 1 and \
                            self(tt[0:len(tt) - 1], trusted=False):
                        yield slice(start, end - 1)
                        pad = min(overlap + 1, end - start - 2)
                        start = end - pad
                    else:
                        # TODO: for some ensembles, there are better ways to
                        # change start. For frame-by-frame ensembles
                        # (AllInX, AllOutX) we know that we can completely
                        # stop for all subtrajectories.
                        start += 1
                    end = start + min_length

        else:
            end = length
            start = end - min_length

            while start >= 0 and end >= min_length:
                tt = trajectory[start:end]

                if len(tt) != old_tt_len + 1:
                    can_prepend_tt = self.can_prepend(tt)
                else:
                    can_prepend_tt = self.can_prepend(tt, trusted=True)
                old_tt_len = len(tt)

                if start > 0 and can_prepend_tt:
                    start -= 1
                    if end - start > max_length + 1:
                        end -= 1
                        start = end - min_length
                else:
                    if end - start <= max_length and self(tt, trusted=False):
                        yield slice(start, end)
                        pad = min(overlap, end - start - 1)
                        end = start + pad
                        if start == 0:
                            # This means we have reached the end and should stop
                            # All other possible subtraj can only be contained
                            # in already existing ones
                            end = 0

                    elif end - start >= min_length + 1 and \
                            self(tt[1:len(tt)], trusted=False):
                        yield slice(start + 1, end)
                        pad = min(overlap + 1, end - start - 2)
                        end = start + pad
                    else:
                        end -= 1

                    start = end - min_length

    def iter_extendable_slices(
            self,
            trajectory,
            max_length=None,
            min_length=1,
            overlap=1,
            reverse=False
    ):
        """
        Return an iterator over maxiaml slices of extendable subtrajectories

        In comparison to the iter_valid_slices this will return maximal
        subtrajectories that can potentially be extended into samples of the
        ensemble. Shorter subparts will also always work. Where we always
        use strict_can_append. So for forward extentable ensembles you can
        cut at the end and for backward extendable ones you can cut at the
        beginning.


        Notes
        -----
        This feature is not yet fully tested and should be used with care!

        Parameters
        ----------
        trajectory : :class:`openpathsampling.trajectory.Trajectory`
            the actual trajectory to be splitted into ensemble parts
        max_length : int > 0, optional
            if set this determines the maximal size to be tested (is mainly
            used in the recursion)
        min_length : int > 0, optional
            if set this determines the minimal size to be tested (in lazy
            mode might no
        overlap : int >= 0, optional
            determines the allowed overlap of all trajectories to be found.
            A value of x means that two sub-trajectorie can share up to x
            frames at the beginning and x frames at the end.  Default is 1
        reverse : bool
            if `True` this will start searching from the end of the trajectory.
            Otherwise (default) it will start at the beginning.

        Returns
        -------
        list of `slice`
            Returns a list of index-slices for sub-trajectories in
            trajectory that are in the ensemble.
        """
        length = len(trajectory)

        logger.info('`iter_extendable_slices` is experimental. Use it on your '
                    'own risk!')

        if max_length is None:
            max_length = length

        max_length = min(length, max_length)
        min_length = max(1, min_length)

        logger.debug("Looking for subtrajectories in " + str(trajectory))
        old_tt_len = 0

        if not reverse:
            start = 0
            end = start + min_length

            while start <= length - min_length and end <= length:
                # print start, end
                tt = trajectory[start:end]

                if len(tt) != old_tt_len + 1:
                    can_append_tt = self.strict_can_append(tt)
                else:
                    can_append_tt = self.strict_can_append(tt, trusted=True)
                old_tt_len = len(tt)

                if end < length and can_append_tt:
                    end += 1
                    if end - start > max_length + 1:
                        start += 1
                        end = start + min_length
                else:
                    if end - start <= max_length + 1:
                        yield slice(start, end - 1)
                        pad = min(overlap, end - start - 1)
                        start = end - pad
                        if end == length:
                            # This means we have reached the end and should stop
                            # All other possible subtraj can only be contained
                            # in already existing ones
                            start = length
                    else:
                        start += 1
                    end = start + min_length

        else:
            end = length
            start = end - min_length

            while start >= 0 and end >= min_length:
                tt = trajectory[start:end]

                if len(tt) != old_tt_len + 1:
                    can_prepend_tt = self.can_prepend(tt)
                else:
                    can_prepend_tt = self.can_prepend(tt, trusted=True)
                old_tt_len = len(tt)

                if start > 0 and can_prepend_tt:
                    start -= 1
                    if end - start > max_length + 1:
                        end -= 1
                        start = end - min_length
                else:
                    if end - start <= max_length + 1:
                        yield slice(start, end - 1)
                        pad = min(overlap, end - start - 1)
                        end = start + pad
                        if start == 0:
                            # This means we have reached the end and should stop
                            # All other possible subtraj can only be contained
                            # in already existing ones
                            end = 0
                    else:
                        end -= 1

                    start = end - min_length

    def find_first_subtrajectory(self, trajectory):
        """
        Return the first sub-trajectory that matches the ensemble

        Parameters
        ----------
        trajectory : :class:`openpathsampling.Trajectory`
            the trajectory in which to look for sub-trajectories

        Returns
        -------
        :class:`openpathsampling.Trajectory` or None
            the found sub-trajectory or None if no sub-trajectory was found
        """
        try:
            return trajectory[
                next(self.iter_valid_slices(trajectory))]
        except StopIteration:
            return None

    def find_last_subtrajectory(self, trajectory):
        """
        Return the last sub-trajectory that matches the ensemble

        Parameters
        ----------
        trajectory : :class:`openpathsampling.trajectory.Trajectory`
            the trajectory in which to look for sub-trajectories

        Returns
        -------
        :class:`openpathsampling.Trajectory` or None
            the found sub-trajectory or None if no sub-trajectory was found
        """
        try:
            return trajectory[
                next(self.iter_valid_slices(trajectory, reverse=True))]
        except StopIteration:
            return None

    def iter_split(
            self,
            trajectory,
            max_length=None,
            min_length=1,
            overlap=1,
            reverse=False):
        """Return iterator over subtrajectories satisfying the given ensemble.

        Parameters
        ----------
        trajectory : :py:class:`openpathsampling.trajectory.Trajectory`
            the actual trajectory to be splitted into ensemble parts
        max_length : int > 0
            if set this determines the maximal size to be tested (is mainly
            used in the recursion)
        min_length : int > 0
            if set this determines the minimal size to be tested (in lazy
            mode might no
        overlap : int >= 0
            determines the allowed overlap of all trajectories to be found.
            A value of x means that two sub-trajectory can share up to x
            frames at the beginning and x frames at the end.  Default is 1
        reverse : bool
            if `True` this will start searching from the end of the trajectory.
            Otherwise (default) it will start at the beginning.

        Returns
        -------
        iterator of :class:`openpathsampling.trajectory.Trajectory`
            Returns a list of sub-trajectories in trajectory that are in the
            ensemble.

        Notes
        -----
        This uses self.iter_valid_slices and returns the actual sub-trajectories
        """
        for part in self.iter_valid_slices(
                trajectory, max_length, min_length, overlap, reverse):
            yield trajectory[part]

    def split(
            self,
            trajectory,
            max_length=None,
            min_length=1,
            overlap=1,
            reverse=False,
            n_results=0):
        """Return list of subtrajectories satisfying the given ensemble.

        Parameters
        ----------
        trajectory : :py:class:`openpathsampling.trajectory.Trajectory`
            the actual trajectory to be splitted into ensemble parts
        max_length : int > 0
            if set this determines the maximal size to be tested (is mainly
            used in the recursion)
        min_length : int > 0
            if set this determines the minimal size to be tested (in lazy
            mode might no
        overlap : int >= 0
            determines the allowed overlap of all trajectories to be found.
            A value of x means that two sub-trajectory can share up to x
            frames at the beginning and x frames at the end.  Default is 1
        reverse : bool
            if `True` this will start searching from the end of the trajectory.
            Otherwise (default) it will start at the beginning.
        n_results : int
            if `0` this will return all results. If the integer is larger than
            zero it will stop after the given number of slices has been found

        Returns
        -------
        list of :class:`openpathsampling.trajectory.Trajectory`
            Returns a list of sub-trajectories in trajectory that are in the
            ensemble.

        Notes
        -----
        This uses self.find_valid_slices and returns the actual sub-trajectories
        """

        indices = self.iter_valid_slices(trajectory, max_length,
                                         min_length, overlap, reverse)

        if n_results > 0:
            return [
                trajectory[part]
                for part in itertools.islice(indices, n_results)]
        else:
            return [trajectory[part] for part in indices]

    @property
    def extendable_sub_ensembles(self):
        return {}

    def get_sample_from_trajectories(
            self, trajectories,
            replica=0,
            used_trajectories=None,
            reuse_strategy='avoid-symmetric'
    ):
        """
        Generate a sample in the ensemble by testing `trajectories`

        Parameters
        ----------
        trajectories : (list of) :class:`openpathsampling.trajectory.Trajectory`
            single trajectory of list of trajectories to be used to create a
            sample in this ensemble
        replica : int
            the replica id for the sample to be created
        used_trajectories : (list of) :class:`openpathsampling.trajectory.Trajectory`
            trajectories not taken into account in the first attempt
        reuse_strategy : str
            if `avoid` then in a second attempt the used trajectories are
            tried
        """

        trajectories = paths.Trajectory._to_list_of_trajectories(trajectories)

        used_and_possible = []

        for idx, traj in enumerate(trajectories):
            if traj not in used_trajectories and (
                    not reuse_strategy.endswith('symmetric') or
                    traj.reversed not in used_trajectories):
                if self(traj):
                    return paths.Sample(
                        trajectory=traj,
                        ensemble=self,
                        replica=replica
                    )
            else:
                used_and_possible.append(traj)

        return self._handle_used_trajectories(
            used_trajectories,
            used_and_possible,
            reuse_strategy)

    def split_sample_from_trajectories(
            self, trajectories,
            replica=0,
            used_trajectories=None,
            reuse_strategy='avoid-symmetric',
            unique='shortest'):
        """
        Generate a sample in the ensemble by searching for sub-parts

        Parameters
        ----------
        trajectories : (list of) :class:`openpathsampling.trajectory.Trajectory`
            single trajectory of list of trajectories to be used to create a
            sample in this ensemble
        replica : int
            the replica id for the sample to be created
        used_trajectories : (list of) :class:`openpathsampling.trajectory.Trajectory`
            trajectories not taken into account in the first attempt
        reuse_strategy : str
            if `avoid` then in a second attempt the used trajectories are
            tried
        unique : str
            If `first` the first found subtrajectory is selected. If
            `shortest` then from all subparts the shortest one is used.
        """

        trajectories = paths.Trajectory._to_list_of_trajectories(trajectories)

        used_and_possible = []

        for idx, traj in enumerate(trajectories):
            parts = self._get_trajectory_parts_in_order(traj, unique)

            for part in parts:
                if part not in used_trajectories and (
                        not reuse_strategy.endswith('symmetric') or
                        part.reversed not in used_trajectories):
                    return paths.Sample(
                        trajectory=part,
                        ensemble=self,
                        replica=replica
                    )
                else:
                    used_and_possible.append(part)

        return self._handle_used_trajectories(
            used_trajectories,
            used_and_possible,
            reuse_strategy)

    def extend_sample_from_trajectories(
            self,
            trajectories,
            engine,
            replica=0,
            unique='median',
            level='complex',
            on_error='retry',
            attempts=2):
        """
        Generate a sample in the ensemble by extending parts of `trajectories`

        This will take an initial trajectory look for useable subparts and
        try to extend them into a valid sample. This works by taking information
        from an ensemble what are resonable subparts, this is returned by a
        function `.extendable_sub_ensembles()` which is only defined for
        complex ensembles like Minus or TIS ensemble.

        As an example the minus could extend from the segment ensemble or even
        a segment + parts completely in the inner ensemble. Of course the
        ensemble itself is always valid.

        The function tries to find extendable subparts from largest to smallest
        ones, starting with the ensemble itself and ending with small subparts

        If a list of trajectories is provided it will be attempt to find a
        valid trajectory using all the trajectory parts.

        Parameters
        ----------
        trajectories : (list of) :class:`openpathsampling.trajectory.Trajectory`
            single trajectory of list of trajectories to be used to create a
            sample in this ensemble
        engine : :class:`openpathsampling.dynamicsengine.DynamicsEngine`
            engine to use for MD extension
        replica : int
            the replica id for the sample to be created
        unique : str
            If `first` the first found subtrajectory is selected. If
            `shortest` then from all subparts the shortest one is used.
        level : str
            there are three levels you chose and not all are implemented for
            an ensemble. For all ensembles you can use `native` which will
            simply try to extend the ensemble itself, the mose simple one, which
            is always possible.
            Picking `complex` will use the largest (most complex)
            sub-ensemble that makes sense. Like in the case of a Minus move
            this is the segment ensemble.
            The other choice is `minimal` which
            choses the minimal necessary subtrajectory extending makes sense
            from. For TIS or Minus Ensembles this will be crossing from the
            (initial) core to the outside. You should try `complex` first and
            then `minimal`. `complex` should be much faster.
        on_error : str
            if `retry` (default) then any error will trigger a retry and
            eventually no sample will be retured. `fail` will raise the
            exception. Typical things to happen are `MaxLengthError` or
            `NaNError`, but also initialisation error can happen. `fail` should
            only be used for debugging purposes since you will not get a
            preliminary sampleset as a result but an exception.
        attempts : int
            the number of attemps on a trajectory to extend
        """

        logger.info("Starting extend_sample_from_trajectories with level "
                    + str(level))
        if level == 'native':
            sub_ensemble = self
        else:
            if not hasattr(self, 'extendable_sub_ensembles'):
                logger.info("Missing ensemble.extendable_sub_ensembles")
                return None

            sub_ensembles = self.extendable_sub_ensembles

            if level not in sub_ensembles:
                logger.info("Missing level: " + repr(level))
                return None

            sub_ensemble = sub_ensembles[level]

        trajectories = paths.Trajectory._to_list_of_trajectories(trajectories)

        for idx, traj in enumerate(trajectories):
            traj_parts = sub_ensemble._get_trajectory_parts_in_order(
                traj, unique)

            for orig in traj_parts:
                for attempt in range(attempts):
                    part = paths.Trajectory(orig)

                    logger.info((
                        'extend - attempt [%d] : extending from initial '
                        'length %d\n') % (
                            attempt + 1,
                            len(part)
                        ))
                    try:
                        if self.strict_can_append(part):
                            # seems we could extend forward

                            part = part[:-1] + \
                               engine.generate(
                                   part[-1],
                                   [paths.PrefixTrajectoryEnsemble(
                                       self,
                                       part
                                   ).strict_can_append],
                                   direction=+1
                               )

                        if self.strict_can_prepend(part):
                            # and extend backward

                            part = engine.generate(
                                part[0].reversed,
                                [paths.SuffixTrajectoryEnsemble(
                                    self,
                                    part
                                ).strict_can_prepend],
                                direction=-1
                            ).reversed + part[1:]

                        logger.info("Candidate trajectory: " + str(part))
                        if self(part):  # make sure we found a sample
                            return paths.Sample(
                                trajectory=part,
                                ensemble=self,
                                replica=replica
                            )
                    except paths.engines.EngineError as e:
                        if on_error == 'fail':
                            raise
                        elif on_error == 'retry':
                            pass
                        else:
                            # This should not happen!
                            pass

        logger.info("Returning None because nothing worked")
        return None

    def _get_trajectory_parts_in_order(self, traj, unique='first'):
        if unique == 'first':
            # this returns an iterator and can thus be faster
            parts = self.iter_split(traj)
        elif unique == 'shortest':
            parts = sorted(self.split(traj), key=len)
        elif unique == 'median':
            # resort the found trajectories so that the middle one is
            # first, then the one right to it, then the one before, etc
            # e.g. [0,1,2,3,4,5,6,7,8,9] is rearranges into
            # [5,4,6,3,7,2,8,1,9,0]
            ordered = sorted(self.split(traj), key=len)
            parts = list([p for p2 in zip(
                ordered[len(ordered) // 2:],
                reversed(ordered[:len(ordered) // 2])
            ) for p in p2])

            if len(ordered) & 1:
                parts.append(ordered[-1])
        elif unique == 'longest':
            parts = sorted(self.split(traj), key=len, reverse=True)
        else:
            parts = []

        try:
            if len(parts) > 0:
                lens = map(len, parts)
                logger.info(
                    ('splitting - found %d slices of lengths '
                     '[%d, ..., %d, ..., %d] '
                     'ordered by `%s`\n') % (
                         len(parts),
                         min(lens),
                         sorted(lens)[len(parts) / 2],
                         max(lens),
                         unique
                     ))
        except TypeError:
            pass

        return parts

    def _handle_used_trajectories(
            self,
            used_trajectories,
            used_and_possible,
            reuse_strategy):

        if reuse_strategy.startswith('avoid') \
                and used_trajectories is not None:

            for part in used_trajectories:
                if part in used_and_possible:
                    if self(part):
                        # move the used one to the back of the list to
                        # not reuse it directly
                        del used_trajectories[used_trajectories.index(part)]
                        used_trajectories.append(part)

                        return paths.Sample(
                            trajectory=part,
                            ensemble=self
                        )

                if reuse_strategy.endswith('symmetric'):
                    if part.reversed in used_and_possible:
                        if self(part):
                            # move the used one to the back of the list to
                            # not reuse it directly
                            del used_trajectories[used_trajectories.index(part)]
                            used_trajectories.append(part)

                            return paths.Sample(
                                trajectory=part,
                                ensemble=self
                            )

        return None

    def __str__(self):
        if self._saved_str is None:
            self._saved_str = self._str()
        return self._saved_str

    def _str(self):
        """
        Returns a complete mathematical expression that defines the current
        ensemble in a readable form.

        Notes
        -----
        This should be cleaned up a little
        """
        return 'Ensemble'

    def __or__(self, other):
        if self is other:
            return self
        elif type(other) is EmptyEnsemble:
            return self
        elif type(other) is FullEnsemble:
            return other
        else:
            return UnionEnsemble(self, other)

            # This is not correct for all ensembles.
            # def __xor__(self, other):
            # # TODO: return (self | other) & ~(self & other)
            # # NOTE: that should also get the automatic special case handling
            # # (other is self, Empty, or Full) from treatment in __and__/__or__
            # if self is other:
            # return EmptyEnsemble()
            # elif type(other) is EmptyEnsemble:
            # return self
            # elif type(other) is FullEnsemble:
            # return NegatedEnsemble(self)
            # else:
            # return SymmetricDifferenceEnsemble(self, other)

    def __and__(self, other):
        if self is other:
            return self
        elif type(other) is EmptyEnsemble:
            return other
        elif type(other) is FullEnsemble:
            return self
        else:
            return IntersectionEnsemble(self, other)

            # This is not correct for all ensembles.
            # def __sub__(self, other):
            # if self is other:
            # return EmptyEnsemble()
            # elif type(other) is EmptyEnsemble:
            # return self
            # elif type(other) is FullEnsemble:
            # return EmptyEnsemble()
            # else:
            # return RelativeComplementEnsemble(self, other)

            # This is not correct for all ensembles.
            # def __invert__(self):
            # return NegatedEnsemble(self)

    @staticmethod
    def _indent(s):
        spl = s.split('\n')
        spl = ['  ' + p for p in spl]
        return '\n'.join(spl)


class EmptyEnsemble(Ensemble):
    """
    The empty path ensemble of no trajectories.
    """

    def __init__(self):
        super(EmptyEnsemble, self).__init__()

    def __call__(self, trajectory, trusted=None, candidate=False):
        return False

    def can_append(self, trajectory, trusted=False):
        return False

    def can_prepend(self, trajectory, trusted=False):
        return False

    def __invert__(self):
        return FullEnsemble()

    def __sub__(self, other):
        return EmptyEnsemble()

    def __and__(self, other):
        return self

    def __xor__(self, other):
        return other

    def __or__(self, other):
        return other

    def _str(self):
        return 'empty'


class FullEnsemble(Ensemble):
    """
    The full path ensemble of all possible trajectories.
    """

    def __init__(self):
        super(FullEnsemble, self).__init__()

    def __call__(self, trajectory, trusted=None, candidate=False):
        return True

    def can_append(self, trajectory, trusted=False):
        return True

    def can_prepend(self, trajectory, trusted=False):
        return True

    def __invert__(self):
        return EmptyEnsemble()

    def __sub__(self, other):
        if type(other) is EmptyEnsemble:
            return self
        elif type(other) is FullEnsemble:
            return EmptyEnsemble()
        else:
            return NegatedEnsemble(other)

    def __and__(self, other):
        return other

    def __xor__(self, other):
        if type(other) is EmptyEnsemble:
            return self
        elif type(other) is FullEnsemble:
            return EmptyEnsemble()
        else:
            return NegatedEnsemble(other)

    def __or__(self, other):
        return self

    def _str(self):
        return 'all'


class NegatedEnsemble(Ensemble):
    """
    Negates an Ensemble and simulates a `not` statement
    """

    # TODO: this whole concept is false and this should be removed
    def __init__(self, ensemble):
        super(NegatedEnsemble, self).__init__()
        self.ensemble = ensemble

    def __call__(self, trajectory, trusted=None, candidate=False):
        return not self.ensemble(trajectory, trusted, candidate)

    def can_append(self, trajectory, trusted=False):
        # We cannot guess the result here so keep on running forever
        return True

    def can_prepend(self, trajectory, trusted=False):
        # We cannot guess the result here so keep on running forever
        return True

    def _str(self):
        return 'not ' + str(self.ensemble)


class EnsembleCombination(Ensemble):
    """
    Logical combination of two ensembles
    """

    def __init__(self, ensemble1, ensemble2, fnc, str_fnc):
        super(EnsembleCombination, self).__init__()
        self.ensemble1 = ensemble1
        self.ensemble2 = ensemble2
        self.fnc = fnc
        self.sfnc = str_fnc
        self.debug = logger.isEnabledFor(logging.DEBUG)

    def to_dict(self):
        return {'ensemble1': self.ensemble1, 'ensemble2': self.ensemble2}

    def _generalized_short_circuit(self, combo, f1, f2, trajectory, trusted,
                                   fname=""):
        """
        Handles short-circuit logic, all in one place for code simplicity.

        Short-circuit logic skips the second part of the combination if the
        result doesn't depend on it.

        Note
        ----
            If you want to enable debug logging for this, it either needs to
            be enabled when the class is instantiated or set with the .debug
            instance variable. This is to improve performance since this
            method is called very frequently.

        Parameters
        ----------
        combo :
            the combination function
        f1 :
            ensemble1's function. Takes trajectory, returns bool. Examples
            include `__call__`, `can_append`, etc.
        f2 :
            ensemble2's function. As with f1, but for ensemble 2.
        trajectory : :class:`.Trajectory`
            input trajectory
        trusted : bool
            the `trusted` flag to send to f1 and f2
        fname : string
            name of the functions f1 and f2. Only used in debug output.
        """
        if self.debug:  # pragma: no cover
            logger.debug("Combination is " + self.__class__.__name__)
        a = f1(trajectory, trusted)
        if self.debug:  # pragma: no cover
            logger.debug("Combination." + fname + ": " +
                         self.ensemble1.__class__.__name__ + " is " + str(a))
            ens2 = f2(trajectory, trusted)
            # logger.debug("Doing ens2_prime")
            # ens2_prime = f2(trajectory, trusted)
            logger.debug("Combination." + fname + ": " +
                         self.ensemble2.__class__.__name__ + " is " + str(ens2))
            # assert(ens2 == ens2_prime)
            logger.debug("Combination should return " + str(self.fnc(a, ens2)))
        res_true = self.fnc(a, True)
        res_false = self.fnc(a, False)
        if res_false == res_true:
            # result is independent of ensemble_b so ignore it
            # logger.debug("Returning res_true == res_false ==" + str(res_true))
            return res_true
        else:
            b = f2(trajectory, trusted)
            # logger.debug("Needs test:" + str(a) + " " + str(self.fnc) +
            #              str(b) + str(self.fnc(a,b)))
            return self.fnc(a, b)

    def __call__(self, trajectory, trusted=None, candidate=False):
        return self._generalized_short_circuit(
            combo=self.fnc,
            f1=self.ensemble1,
            f2=self.ensemble2,
            trajectory=trajectory,
            trusted=trusted,
            fname="__call__"
        )

    def can_append(self, trajectory, trusted=False):
        return self._generalized_short_circuit(
            combo=self.fnc,
            f1=self.ensemble1.can_append,
            f2=self.ensemble2.can_append,
            trajectory=trajectory,
            trusted=trusted,
            fname="can_append"
        )

    def can_prepend(self, trajectory, trusted=False):
        return self._generalized_short_circuit(
            combo=self.fnc,
            f1=self.ensemble1.can_prepend,
            f2=self.ensemble2.can_prepend,
            trajectory=trajectory,
            trusted=trusted,
            fname="can_prepend"
        )

    def strict_can_append(self, trajectory, trusted=False):
        return self._generalized_short_circuit(
            combo=self.fnc,
            f1=self.ensemble1.strict_can_append,
            f2=self.ensemble2.strict_can_append,
            trajectory=trajectory,
            trusted=trusted,
            fname="strict_can_append"
        )

    def strict_can_prepend(self, trajectory, trusted=False):
        return self._generalized_short_circuit(
            combo=self.fnc,
            f1=self.ensemble1.strict_can_prepend,
            f2=self.ensemble2.strict_can_prepend,
            trajectory=trajectory,
            trusted=trusted,
            fname="strict_can_prepend"
        )

    def _str(self):
        # print self.sfnc, self.ensemble1, self.ensemble2,
        # print self.sfnc.format(
        #     '(' + str(self.ensemble1) + ')',
        #     '(' + str(self.ensemble1) + ')')
        return self.sfnc.format(
            '(\n' + Ensemble._indent(str(self.ensemble1)) + '\n)',
            '(\n' + Ensemble._indent(str(self.ensemble2)) + '\n)')


class UnionEnsemble(EnsembleCombination):
    def __init__(self, ensemble1, ensemble2):
        super(UnionEnsemble, self).__init__(ensemble1, ensemble2,
                                            fnc=lambda a, b: a or b,
                                            str_fnc='{0}\nor\n{1}')


class IntersectionEnsemble(EnsembleCombination):
    def __init__(self, ensemble1, ensemble2):
        super(IntersectionEnsemble, self).__init__(ensemble1, ensemble2,
                                                   fnc=lambda a, b: a and b,
                                                   str_fnc='{0}\nand\n{1}')


# class SymmetricDifferenceEnsemble(EnsembleCombination):
#     # TODO: this is not yet supported. Should be removed. ~DWHS
#     # should just be a shortcut for (ens1 | ens2) & ~(ens1 & ens2)
#     # should probably not even be a class. Just have `ensemble.__xor__`
#     # return (ens1 | ens2) & ~(ens1 & ens2)
#     def __init__(self, ensemble1, ensemble2):
#         super(SymmetricDifferenceEnsemble, self).__init__(
#             ensemble1,
#             ensemble2,
#             fnc=lambda a, b: a ^ b,
#             str_fnc='{0}\nxor\n{1}')


# class RelativeComplementEnsemble(EnsembleCombination):
#     # TODO: this is not yet supported. Should be removed. ~DWHS
#     # should be a shortcut for ens1 & ~ens2
#     # should probably not even be a class. Just have `ensemble.__sub__`
#     # return ens1 & ~ens2
#     def __init__(self, ensemble1, ensemble2):
#         super(RelativeComplementEnsemble, self).__init__(
#             ensemble1,
#             ensemble2,
#             fnc=lambda a, b: a and not b,
#             str_fnc='{0}\nand not\n{1}')


class SequentialEnsemble(Ensemble):
    """Ensemble which satisfies several subensembles in sequence.

    Attributes
    ----------
    ensembles : tuple of Ensemble
        The ensembles, in time-order of when they should occur in the
        trajectory.
    min_overlap : int or tuple of int
        The minimum number of frames that overlap between two ensembles in
        the sequence. A positive number n indicates that at least n frames
        must be in both ensembles at the transition between them. A negative
        number -n indicates that at least n frames in neither ensemble at
        the transition between them. If given as a list, the list should be
        of length len(ensembles)-1, with one value for each transition. If
        given as an integer, that value will be used for all transitions.
    max_overlap : int or list of int
        The maximum number of frames that overlap between two ensembles in
        the sequence. A positive number n indicates that no more than n
        frames can be in both ensembles at the transition between them. A
        negative number -n indicates no more than n frames in neither
        ensemble at the transition between them. If given as a list, the
        list should be of length len(ensembles)-1, with one value for each
        transition. If given as an integer, that value will be used for all
        transitions.

    Notes
    -----
        TODO: Overlap features not implemented because ohmygod this was hard
        enough already.
    """

    def __init__(self, ensembles, min_overlap=0, max_overlap=0, greedy=False):
        # make tuples of the min/max overlaps
        super(SequentialEnsemble, self).__init__()
        if type(min_overlap) is int:
            min_overlap = (min_overlap,) * (len(ensembles) - 1)
        if type(max_overlap) is int:
            max_overlap = (max_overlap,) * (len(ensembles) - 1)

        self.ensembles = ensembles
        self.min_overlap = min_overlap
        self.max_overlap = max_overlap
        self.greedy = greedy

        self._use_cache = True  # cache can be turned off
        self._cache_can_append = EnsembleCache(+1)
        self._cache_strict_can_append = EnsembleCache(+1)
        self._cache_call = EnsembleCache(+1)
        self._cache_can_prepend = EnsembleCache(-1)
        self._cache_strict_can_prepend = EnsembleCache(-1)
        self._cache_check_reverse = EnsembleCache(-1)
        self._zero_traj = paths.Trajectory([])

        # sanity checks
        if len(self.min_overlap) != len(self.max_overlap):
            raise ValueError("len(min_overlap) != len(max_overlap)")
        if len(self.min_overlap) != len(self.ensembles) - 1:
            raise ValueError(
                "Number of overlaps doesn't match number of transitions")
        for i in range(len(self.min_overlap)):
            if min_overlap[i] > max_overlap[i]:
                raise ValueError("min_overlap greater than max_overlap!")

    @staticmethod
    def update_cache(cache, ens_num, ens_from, subtraj_from):
        """Updates the given cache.

        Parameters
        ----------
        cache : `EnsembleCache`
            the cache to be updated
        ens_num : integer
            current value of `ens_num` in the sequential ensemble
        ens_from : integer
            current "start" ensemble index. For forward-direction caches,
            this is ens_first. For reverse-direction caches, this is
            ens_final. The "initial" (in the appropriate direction) frame is
            assigned to this ensemble
        subtraj_from : integer
            index of the "start" frame of the subtrajectory in this
            subensemble. For forward-direction caches, this is the first
            frame of the subtrajectory. For reverse-direction caches, this
            is the final frame of the subtrajectory.
        """
        if ens_num == "keep":
            ens_num = cache.contents['ens_num']
        if ens_from == "keep":
            ens_from = cache.contents['ens_from']
        if subtraj_from == "keep":
            subtraj_from = cache.contents['subtraj_from']

        cache.contents['ens_num'] = ens_num
        cache.contents['ens_from'] = ens_from
        cache.contents['subtraj_from'] = subtraj_from
        logger.debug("Setting cache | ens_num " + str(ens_num) +
                     " | ens_from " + str(ens_from) +
                     " | subtraj_from " + str(subtraj_from))
        logger.debug("Cache is Trusted: " + str(cache.trusted))

    @staticmethod
    def assign_frames(cache, ens_num,
                      subtraj_first=None, subtraj_final=None):
        if ens_num is None:
            cache.contents['assignments'] = {}
        else:
            cache.contents['assignments'][ens_num] = \
                slice(subtraj_first, subtraj_final)
        logger.debug("Cache assignments: " + str(cache.contents['assignments']))

    def transition_frames(self, trajectory, trusted=None):
        # it is easiest to understand this decision tree as a simplified
        # version of the can_append decision tree; see that for detailed
        # comments
        # self._check_cache(trajectory, function="call")

        ens_num = 0
        subtraj_first = 0

        traj_final = len(trajectory)
        final_ens = len(self.ensembles) - 1
        transitions = []
        while True:
            if ens_num <= final_ens:
                subtraj_final = self._find_subtraj_final(trajectory,
                                                         subtraj_first, ens_num)
            else:
                return transitions
            if subtraj_final - subtraj_first > 0:
                # subtraj = trajectory[slice(subtraj_first, subtraj_final)]
                if ens_num == final_ens:
                    if subtraj_final == traj_final:
                        # success
                        transitions.append(subtraj_final)
                        return transitions
                    else:
                        # fails because we have more frames to assign
                        transitions.append(subtraj_final)
                        return transitions
                else:
                    ens_num += 1
                    transitions.append(subtraj_final)
                    subtraj_first = subtraj_final
            else:
                if ens_num <= final_ens and \
                        self.ensembles[ens_num](self._zero_traj):
                    ens_num += 1
                    transitions.append(subtraj_final)
                    subtraj_first = subtraj_final
                else:
                    return transitions

    def __call__(self, trajectory, trusted=None, candidate=False):
        logger.debug("Looking for transitions in trajectory " + str(trajectory))
        transitions = self.transition_frames(trajectory, trusted)
        logger.debug("Found transitions: " + str(transitions))
        # if we don't have the right number of transitions, or if the last
        # print transitions
        if len(transitions) != len(self.ensembles):
            # print "Returns false b/c not enough ensembles"
            return False
        elif transitions[-1] != len(trajectory):
            # print "Returns false b/c not all frames assigned"
            return False

        subtraj_first = 0
        subtraj_i = 0
        # Make a list before slicing
        ltraj = _get_list_traj(trajectory)
        while subtraj_i < len(self.ensembles):
            subtraj_final = transitions[subtraj_i]
            subtraj = ltraj[slice(subtraj_first, subtraj_final)]
            if not self.ensembles[subtraj_i](subtraj):
                # print "Returns false b/c ensemble", subtraj_i," fails"
                return False
            subtraj_i += 1
            subtraj_first = subtraj_final
        return True

    def _find_subtraj_final(self, traj, subtraj_first, ens_num,
                            last_checked=None):
        """
        Find the longest subtrajectory of trajectory which starts at
        subtraj_first and satifies self.ensembles[ens_num].can_append

        Returns
        -------
        int
            Frame of traj which is the final frame for a subtraj starting at
            subtraj_first and satisfying self.ensembles.can_append[ens_num]
        """
        if last_checked is None:
            subtraj_final = subtraj_first
        else:
            subtraj_final = max(last_checked, subtraj_first)
        traj_final = len(traj)
        ens = self.ensembles[ens_num]
        # Make a list before slicing
        ltraj = _get_list_traj(traj)
        subtraj = ltraj[slice(subtraj_first, subtraj_final + 1)]

        # if we're in the ensemble or could eventually be in the ensemble,
        # we keep building the subtrajectory

        # TODO: this doesn't actually reflect the cleanest behavior: should
        # be the proper hybrid definition where we can append until/unless
        # we overshoot
        logger.debug("*Traj slice " + str(subtraj_first) + " " +
                     str(subtraj_final + 1) + " / " + str(traj_final))
        # logger.debug("Ensemble " + str(ens.__class__.__name__))# + str(ens))
        # logger.debug("Can-app " + str(ens.can_append(subtraj, trusted=True)))
        # logger.debug("Call    " + str(ens(subtraj, trusted=True)))
        # TODO: the weird while condition is handling the OVERSHOOTING
        while ((ens.can_append(subtraj, trusted=True) or
                ens(subtraj, trusted=True)) and subtraj_final < traj_final):
            subtraj_final += 1
            subtraj = ltraj[slice(subtraj_first, subtraj_final + 1)]
            logger.debug(" Traj slice " + str(subtraj_first) + " " +
                         str(subtraj_final + 1) + " / " + str(traj_final))
        return subtraj_final

    def _find_subtraj_first(self, traj, subtraj_final, ens_num,
                            last_checked=None):
        if last_checked is None:
            subtraj_first = subtraj_final - 1
        else:
            subtraj_first = min(last_checked, subtraj_final - 1)
        traj_first = 0
        ens = self.ensembles[ens_num]

        # Make a list before slicing
        ltraj = _get_list_traj(traj)
        subtraj = ltraj[slice(subtraj_first, subtraj_final)]
        logger.debug("*Traj slice " + str(subtraj_first) + " " +
                     str(subtraj_final) + " / " + str(len(traj)))
        # logger.debug("Ensemble " + str(ens.__class__.__name__))# + str(ens))
        # logger.debug("Can-app " + str(ens.can_prepend(subtraj, trusted=True)))
        # logger.debug("Call    " + str(ens(subtraj, trusted=True)))
        # TODO: the weird while condition is handling the OVERSHOOTING
        while ((ens.can_prepend(subtraj, trusted=True) or
                ens.check_reverse(subtraj, trusted=True)
               ) and subtraj_first >= traj_first):
            subtraj_first -= 1
            subtraj = ltraj[slice(subtraj_first, subtraj_final)]
            logger.debug(" Traj slice " + str(subtraj_first + 1) + " " +
                         str(subtraj_final) + " / " + str(len(traj)))
        return subtraj_first + 1

    def _generic_can_append(self, trajectory, trusted, strict):
        # treat this like we're implementing a regular expression parser ...
        # .*ensemble.+ ; but we have to do this for all possible matches
        # There are three tests we consider:
        # 1. subtraj_final - subtraj_first > 0: Do we obtain a subtrajectory?
        # 2. subtraj_final == traj_final: Have we assigned all the frames?
        # 3. ens_num == final_ens: are we looking at the last ensemble
        # Various combinations of these result in three possible outcomes:
        # (a) return True (we can append)
        # (b) return False (we can't append)
        # (c) loop around to text another subtrajectory (we can't tell)
        # Returning false can only happen if all ensembles have been tested
        # self._check_cache(trajectory, function="can_append")
        cache = self._cache_can_append
        if strict:
            cache = self._cache_strict_can_append

        if trusted:
            cache.trusted = True

        subtraj_first = 0
        ens_num = 0
        ens_first = 0

        if self._use_cache:
            _ = cache.check(trajectory)
            if cache.contents == {}:
                self.update_cache(cache, 0, 0, 0)
                self.assign_frames(cache, None)
            else:
                subtraj_first = cache.contents['subtraj_from']
                ens_num = cache.contents['ens_num']
                ens_first = cache.contents['ens_from']

        traj_final = len(trajectory)
        final_ens = len(self.ensembles) - 1
        # Make a list before slicing
        ltraj = _get_list_traj(trajectory)
        # print traj_final, final_ens
        # logging startup
        if cache.debug_enabled:  # pragma: no cover
            logger.debug(
                "Beginning can_append with subtraj_first="
                + str(subtraj_first) + "; ens_first=" + str(ens_first)
                + "; ens_num=" + str(ens_num)
                + "; strict=" + str(strict)
            )
            logger.debug(
                "Can-append sees a trusted cache: " + str(cache.trusted)
            )
            if cache.trusted:
                logger.debug("Cache contents: " + str(cache.contents))
                logger.debug("cache.prev_last_frame: " +
                             str(trajectory.index(cache.prev_last_frame)))
            for i in range(len(self.ensembles)):
                ens = self.ensembles[i]
                logger.debug("Ensemble " + str(i) + " : "
                             + ens.__class__.__name__)

        while True:  # main loop, with various exits
            if self._use_cache and cache.trusted:
                # TODO: trajectory.index is expensive... how to speed up?
                # offset = 1
                offset = 0
                # if cache.last_length == len(trajectory):
                # offset += 1
                try:
                    last_checked_index = cache.prev_last_index - offset
                except:  # on any exception
                    # TODO: ideally, this won't be covered by tests, and can
                    # eventually be removed (along with try/except)
                    last_checked_index = \
                        trajectory.index(cache.prev_last_frame) - offset
                #last_checked = trajectory.index(cache.prev_last_frame) - offset
                last_checked = last_checked_index
            else:
                last_checked_index = None
                last_checked = None
            if cache.debug_enabled:
                logger.debug("last_checked = " + str(last_checked))
            subtraj_final = self._find_subtraj_final(
                trajectory, subtraj_first, ens_num, last_checked
            )
            cache.last_length = subtraj_final
            if cache.debug_enabled:
                logger.debug(
                    "Subtraj for ens " + str(ens_num) + " : " +
                    "(" + str(subtraj_first) + "," + str(subtraj_final) + ")"
                )
            if subtraj_final - subtraj_first > 0:
                subtraj = ltraj[slice(subtraj_first, subtraj_final)]
                if ens_num == final_ens:
                    if subtraj_final == traj_final:
                        # we're in the last ensemble and the whole
                        # trajectory is assigned: can we append?
                        ens = self.ensembles[ens_num]
                        if cache.debug_enabled:
                            logger.debug("Returning can_append for " +
                                         str(ens.__class__.__name__))
                        self.update_cache(cache, ens_num,
                                          ens_first, subtraj_first)
                        return ens.can_append(subtraj, trusted=True)
                    else:
                        logger.debug(
                            "Returning false due to incomplete assigns: " +
                            str(subtraj_final) + "!=" + str(traj_final)
                        )
                        return False  # in final ensemble, not all assigned
                else:
                    # subtraj existed, but not yet final ensemble
                    # so we start with the next ensemble
                    end_traj = (subtraj_final == traj_final)
                    ensemble = self.ensembles[ens_num]
                    if not end_traj and not ensemble(subtraj,
                                                     trusted=cache.trusted):
                        logger.debug(
                            "Couldn't assign frames " + str(subtraj_first) +
                            " through " + str(subtraj_final) +
                            " to ensemble " + str(ens_num) + ": No match"
                        )
                    else:
                        logger.debug(
                            "Assigning frames " + str(subtraj_first) +
                            " through " + str(subtraj_final) +
                            " to ensemble " + str(ens_num)
                        )
                        self.assign_frames(cache, ens_num, subtraj_first,
                                           subtraj_final)
                        self.update_cache(cache, ens_num, ens_first,
                                          subtraj_first)
                    ens_num += 1
                    subtraj_first = subtraj_final
                    logger.debug("Moving to the next ensemble " + str(ens_num))
            else:  # no subtrajectory found
                if subtraj_final == traj_final:
                    # all frames assigned, but not all ensembles finished;
                    # next frame might satisfy next ensemble
                    if self._use_cache:
                        prev_slice = cache.contents['assignments'][ens_num - 1]
                        prev_subtraj = ltraj[prev_slice]
                        prev_ens = self.ensembles[ens_num - 1]
                        if prev_ens.can_append(prev_subtraj, trusted=True):
                            logger.debug(
                                "Premature promotion: returning to ensemble " +
                                str(ens_num - 1)
                            )
                            ens_num -= 1
                            subtraj_first = "keep"

                        self.update_cache(cache, ens_num, ens_first,
                                          subtraj_first)
                    logger.debug(
                        "All frames assigned, more ensembles to go: "
                        "returning True")
                    return True

                elif self.ensembles[ens_num](self._zero_traj):
                    logger.debug(
                        "Moving on because of allowed zero-length ensemble")
                    ens_num += 1
                    subtraj_first = subtraj_final
                    self.update_cache(cache, ens_num, ens_first, subtraj_first)

                else:
                    # not all frames assigned, couldn't find a sequence
                    # start over with sequences that begin with the next
                    # ensemble
                    if ens_first == final_ens:
                        logger.debug(
                            "Started with the last ensemble, got nothin'")
                        return False
                    elif strict is False:
                        logger.debug(
                            "Reassigning all frames, starting with ensemble " +
                            str(ens_first)
                        )
                        ens_first += 1
                        ens_num = ens_first
                        subtraj_first = 0
                        self.update_cache(cache, ens_num, ens_first,
                                          subtraj_first)
                    else:
                        logger.debug(
                            "First ensemble fails and strict -- return false"
                        )
                        return False

    def can_append(self, trajectory, trusted=False):
        return self._generic_can_append(trajectory, trusted, strict=False)

    def strict_can_append(self, trajectory, trusted=False):
        return self._generic_can_append(trajectory, trusted, strict=True)

    def _generic_can_prepend(self, trajectory, trusted, strict):
        # based on .can_append(); see notes there for algorithm details
        cache = self._cache_can_prepend
        if strict:
            cache = self._cache_strict_can_prepend
        if trusted:
            cache.trusted = True

        traj_first = 0
        first_ens = 0
        subtraj_final = len(trajectory)
        ens_final = len(self.ensembles) - 1
        ens_num = ens_final
        # Make list before slicing
        ltraj = _get_list_traj(trajectory)

        if self._use_cache:
            _ = cache.check(trajectory)
            if cache.contents == {}:
                self.update_cache(cache, ens_num, first_ens, subtraj_final)
                self.assign_frames(cache, None)
            else:
                logger.debug(
                    "len(traj)=" + str(len(trajectory))
                    + "cache_from=" + str(cache.contents['subtraj_from'])
                )
                subtraj_from = cache.contents['subtraj_from']
                if subtraj_from is None:
                    subtraj_from = 0
                subtraj_final = len(trajectory) + subtraj_from
                ens_num = cache.contents['ens_num']
                ens_final = cache.contents['ens_from']

        # logging startup
        if logger.isEnabledFor(logging.DEBUG):  # pragma: no cover
            logger.debug(
                "Beginning can_prepend with ens_num:" + str(ens_num)
                + "  ens_final:" + str(ens_final) + "  subtraj_final "
                + str(subtraj_final) + "; strict=" + str(strict)
            )
            if cache.trusted:
                logger.debug("Cache contents: " + str(cache.contents))
                logger.debug("cache.prev_start_frame: " +
                             str(trajectory.index(cache.start_frame)))
            for i in range(len(self.ensembles)):
                logger.debug(
                    "Ensemble " + str(i) +
                    " : " + self.ensembles[i].__class__.__name__
                )

        while True:
            if self._use_cache and cache.trusted:
                # offset = 1
                offset = 0
                try:
                    last_checked_index = cache.prev_last_index + offset
                except:  # on any exception
                    # TODO: ideally, this won't be covered by tests, and can
                    # eventually be removed (along with try/except)
                    last_checked_index = \
                        trajectory.index(cache.prev_last_frame) + offset
                last_checked = trajectory.index(cache.prev_last_frame) + offset
            else:
                last_checked_index = None
                last_checked = None
            subtraj_first = self._find_subtraj_first(
                trajectory, subtraj_final, ens_num, last_checked)
            cache.last_length = len(trajectory) - subtraj_first

            assign_final = subtraj_final - len(trajectory)
            if assign_final == 0:
                assign_final = None
            logger.debug(
                str(ens_num) + " : " +
                "(" + str(subtraj_first) + "," + str(subtraj_final) + ")"
            )
            if subtraj_final - subtraj_first > 0:
                subtraj = ltraj[slice(subtraj_first, subtraj_final)]
                if ens_num == first_ens:
                    if subtraj_first == traj_first:
                        logger.debug("Returning can_prepend")
                        self.update_cache(cache, ens_num, ens_final,
                                          assign_final)
                        return self.ensembles[ens_num].can_prepend(subtraj,
                                                                   trusted=True)
                    else:
                        logger.debug(
                            "Returning false due to incomplete assigns: " +
                            str(subtraj_first) + "!=" + str(traj_first)
                        )
                        return False
                else:
                    if subtraj_first != traj_first and \
                            not self.ensembles[ens_num](
                                subtraj, trusted=True):
                        logger.debug(
                            "Couldn't assign frames " + str(subtraj_first) +
                            " through " + str(subtraj_final) +
                            " to ensemble " + str(ens_num) + ": No match"
                        )
                    else:
                        logger.debug(
                            "Assigning frames " + str(subtraj_first) +
                            " through " + str(subtraj_final) +
                            " to ensemble " + str(ens_num)
                        )
                        assign_first = subtraj_first - len(trajectory)
                        self.assign_frames(cache, ens_num, assign_first,
                                           assign_final)
                        self.update_cache(cache, ens_num, ens_final,
                                          assign_final)
                    ens_num -= 1
                    subtraj_final = subtraj_first
                    logger.debug("Moving to the next ensemble " + str(ens_num))
            else:
                if subtraj_first == traj_first:
                    if self._use_cache:
                        prev_slice = cache.contents['assignments'][ens_num + 1]
                        logger.debug("prev_slice " + str(prev_slice))
                        prev_subtraj = ltraj[prev_slice]
                        logger.debug("prev_subtraj " + str(prev_subtraj))
                        logger.debug("traj " + str(trajectory))
                        prev_ens = self.ensembles[ens_num + 1]
                        if prev_ens.can_prepend(prev_subtraj, trusted=True):
                            logger.debug(
                                "Premature promotion: returning to ensemble " +
                                str(ens_num + 1)
                            )
                            ens_num += 1
                            assign_final = "keep"

                        logger.debug("(first, final)" + str((subtraj_first,
                                                             subtraj_final)))
                        self.update_cache(cache, ens_num, ens_final,
                                          assign_final)
                    logger.debug(
                        "All frames assigned, more ensembles to go: "
                        "returning True")
                    return True
                elif self.ensembles[ens_num](self._zero_traj):
                    logger.debug(
                        "Moving on because of allowed zero-length ensemble")
                    ens_num -= 1
                    subtraj_final = subtraj_first
                    self.update_cache(cache, ens_num, ens_final, subtraj_final)
                else:
                    if ens_final == first_ens:
                        logger.debug(
                            "Started with the last ensemble, got nothin'")
                        return False
                    elif strict is False:
                        logger.debug(
                            "Reassigning all frames, starting with ensemble " +
                            str(ens_final)
                        )
                        ens_final -= 1
                        ens_num = ens_final
                        subtraj_final = len(trajectory)
                        self.update_cache(cache, ens_num, ens_final,
                                          subtraj_final)
                    else:
                        logger.debug(
                            "First ensemble fails and strict -- return false"
                        )
                        return False

    def can_prepend(self, trajectory, trusted=False):
        return self._generic_can_prepend(trajectory, trusted, strict=False)

    def strict_can_prepend(self, trajectory, trusted=False):
        return self._generic_can_prepend(trajectory, trusted, strict=True)

    def _str(self):
        head = "[\n"
        tail = "\n]"
        sequence_str = ",\n".join([str(ens) for ens in self.ensembles])
        return head + sequence_str + tail


class LengthEnsemble(Ensemble):
    """
    The ensemble of trajectories of a given length
    """

    def __init__(self, length):
        """
        A path ensemble that describes path of a specific length

        Parameters
        ----------
        length : int or slice
            The specific length (int) or the range of allowed trajectory
            lengths (slice)
        """
        #TODO: remove support for slice?

        super(LengthEnsemble, self).__init__()
        self.length = length

    def __call__(self, trajectory, trusted=None, candidate=False):
        length = len(trajectory)
        if type(self.length) is int:
            return length == self.length
        else:
            return length >= self.length.start and (
                self.length.stop is None or length < self.length.stop)

    def can_append(self, trajectory, trusted=False):
        length = len(trajectory)
        if type(self.length) is int:
            return_value = (length < self.length)
            logger.debug("LengthEnsemble.can_append: Segment length " +
                         str(length) + " < " + str(self.length) + " : " +
                         str(return_value))
            return return_value
        else:
            return self.length.stop is None or length < self.length.stop - 1

    def can_prepend(self, trajectory, trusted=False):
        return self.can_append(trajectory)

    def _str(self):
        if type(self.length) is int:
            return 'len(x) = {0}'.format(self.length)
        else:
            start = self.length.start
            if start is None:
                start = 0
            stop = self.length.stop
            if stop is None:
                stop = 'infty'
            else:
                stop = str(self.length.stop - 1)
            return 'len(x) in [{0}, {1}]'.format(start, stop)


class VolumeEnsemble(Ensemble):
    """
    Path ensembles based on the Volume object
    """

    def __init__(self, volume, trusted=True):
        # TODO: does `trusted` actually mean anything or do anything as a
        # property? it is about the condition of trusting the trajectory
        # when we run it, so it relevant in functions. I don't think we need
        # it here. ~DWHS
        super(VolumeEnsemble, self).__init__()
        self.volume = volume
        self.trusted = trusted

        self._use_cache = True
        self._cache_can_append = EnsembleCache(+1)
        self._cache_call = EnsembleCache(+1)
        self._cache_can_prepend = EnsembleCache(-1)
        self._cache_check_reverse = EnsembleCache(-1)

    @property
    def _volume(self):
        """
        The volume that is used in the specification.
        """
        return self.volume


class AllInXEnsemble(VolumeEnsemble):
    """
    Ensemble of trajectories with all frames in the given volume
    """

    def _trusted_call(self, trajectory, cache):
        """
        Generalized version of the call when trusted.

        This uses a cache, which has the result for the previous trajectory
        (`trajectory[:-1]` if forward, `trajectory[1:]` if backward) in the
        `cache.contents['previous']`.

        Paramters
        ---------
        trajectory : paths.Trajectory
            input trajectory to test
        cache : paths.EnsembleCache
            ensemble cache for this function

        Returns
        -------
        bool :
            result of __call__
        """
        frame_num = -(cache.direction + 1) // 2  # 1 -> -1; -1 -> 0
        reset = cache.check(trajectory)
        if reset:
            if len(trajectory) < 2:
                cache.contents['previous'] = None
            else:
                # NOTE: is it possible that we'd reset a cache more than
                # once in a single trajectory? that could mean that this
                # starts to scale quadratically. I can't think of a case
                # where this is a practical concern (short-circuit logic
                # means the recache should only happen once per trajectory
                # for All*XEnsembles, and the call should only happen once
                # per trajectory for Part*XEnsembles.) In any case, the fix
                # would be to implement a more complicated cache.reset,
                # which checks whether the previous traj was a subtraj of
                # this one (other than one frame less). ~~~DWHS
                if frame_num == -1:
                    reset_value = self(trajectory[:-1], trusted=False)
                elif frame_num == 0:
                    reset_value = self(trajectory[1:], trusted=False)
                else:  # pragma: no cover
                    raise RuntimeError("Bad value for frame_num: " +
                                       str(frame_num))
                cache.contents['previous'] = reset_value

        cached_val = cache.contents['previous']
        if cached_val or cached_val is None:
            # need to check this frame (no prev traj, or prev traj is True)
            # This sometimes gets a list instead of a full Trajectory
            get_frame = getattr(trajectory, "get_as_proxy",
                                trajectory.__getitem__)
            frame = get_frame(frame_num)
            cache.contents['previous'] = self._volume(frame)
            return cache.contents['previous']
        else:
            # cached_val is false, result must be false
            return False

    def can_append(self, trajectory, trusted=False):
        if len(trajectory) == 0:
            return True
        elif trusted and self._use_cache:
            return self._trusted_call(trajectory, self._cache_can_append)
        else:
            return self(trajectory)

    def can_prepend(self, trajectory, trusted=False):
        if len(trajectory) == 0:
            return True
        if trusted and self._use_cache:
            return self._trusted_call(trajectory, self._cache_can_prepend)
        else:
            return self(trajectory)

    def __call__(self, trajectory, trusted=None, candidate=False):
        if len(trajectory) == 0:
            return False
        # TODO: We might be able to speed this up based on can_append
        # being the same as call for this ensemble. Something like check
        # the can_append cache instead of/as well as the call cache. May
        # still have problems with overshooting -- but this might provide a
        # speed-up in sequential ensemble's checking phase. ~~~DWHS
        if trusted and self._use_cache:
            return self._trusted_call(trajectory, self._cache_call)
        else:
            logger.debug("Untrusted VolumeEnsemble " + repr(self))
            # logger.debug("Trajectory " + repr(trajectory))
            # This can sometimes get a list instead of a Trajectory
            # Make sure this is a proxy list
            for frame in _get_list_traj(trajectory):
                if not self._volume(frame):
                    return False
            return True

    def check_reverse(self, trajectory, trusted=False):
        # order in this one only matters if it is trusted
        if trusted and self._use_cache:
            # print "Rev Trusted"
            return self._trusted_call(trajectory, self._cache_check_reverse)
            # frame = trajectory.get_as_proxy(0)
            # return self._volume(frame)
        else:
            # print "Rev UnTrusted"
            return self(trajectory)  # in this case, order wouldn't matter

    def __invert__(self):
        return PartOutXEnsemble(self.volume, self.trusted)

    def _str(self):
        return 'x[t] in {0} for all t'.format(self._volume)


class AllOutXEnsemble(AllInXEnsemble):
    """
    Ensemble of trajectories with all frames outside the given volume
    """
    def __init__(self, volume, trusted=True):
        super(AllOutXEnsemble, self).__init__(volume, trusted)
        self._cached_volume = ~self.volume

    @property
    def _volume(self):
        return self._cached_volume

    def _str(self):
        return 'x[t] in {0} for all t'.format(self._volume)

    def __invert__(self):
        return PartInXEnsemble(self.volume, self.trusted)


class PartInXEnsemble(VolumeEnsemble):
    """
    Ensemble of trajectory with at least one frame in the volume
    """

    def _str(self):
        return 'exists t such that x[t] in {0}'.format(self._volume)

    def __call__(self, trajectory, trusted=None, candidate=False):
        """
        Returns True if the trajectory is part of the PathEnsemble

        Parameters
        ----------
        trajectory : :class:`openpathsampling.trajectory.Trajectory`
            The trajectory to be checked
        """
        for frame in _get_list_traj(trajectory):
            if self._volume(frame):
                return True
        return False

    def __invert__(self):
        return AllOutXEnsemble(self.volume, self.trusted)


class PartOutXEnsemble(PartInXEnsemble):
    """
    Ensemble of trajectories with at least one frame outside the volume
    """
    def __init__(self, volume, trusted=True):
        super(PartOutXEnsemble, self).__init__(volume, trusted)
        self._cached_volume = ~self.volume

    @property
    def _volume(self):
        return self._cached_volume

    def _str(self):
        return 'exists t such that x[t] in {0}'.format(self._volume)

    def __invert__(self):
        return AllInXEnsemble(self.volume, self.trusted)

    def __call__(self, trajectory, trusted=None, candidate=False):
        # Don't load proxies if this is a Trajectory
        for frame in _get_list_traj(trajectory):
            if self._volume(frame):
                return True
        return False


class WrappedEnsemble(Ensemble):
    """
    Wraps an ensemble to alter it or the way it sees a trajectory
    """

    def __init__(self, ensemble):
        super(WrappedEnsemble, self).__init__()
        self.ensemble = ensemble

        # you can also build wrapped ensembles with more flexibility when using
        # a property for _new_ensemble
        self._new_ensemble = self.ensemble
        self.trusted = None
        self._cache_can_append = EnsembleCache(+1)
        self._cache_strict_can_append = EnsembleCache(+1)
        self._cache_call = EnsembleCache(+1)

        # cache_can_prepend has to think it is going forward because the
        # frames given to it are from a forward growing trajectory... only
        # later is everything turned around
        self._cache_can_prepend = EnsembleCache(+1)
        self._cache_strict_can_prepend = EnsembleCache(+1)

    def __call__(self, trajectory, trusted=None, candidate=False):
        return self._new_ensemble(self._alter(trajectory), trusted)

    def _alter(self, trajectory):
        return trajectory

    def can_append(self, trajectory, trusted=None):
        return self._new_ensemble.can_append(self._alter(trajectory),
                                             trusted)

    def can_prepend(self, trajectory, trusted=None):
        return self._new_ensemble.can_prepend(self._alter(trajectory),
                                              trusted)

    def strict_can_append(self, trajectory, trusted=None):
        return self._new_ensemble.strict_can_append(self._alter(trajectory),
                                                    trusted)

    def strict_can_prepend(self, trajectory, trusted=None):
        return self._new_ensemble.strict_can_prepend(self._alter(trajectory),
                                                     trusted)

    def _str(self):
        return str(self._new_ensemble)


class SlicedTrajectoryEnsemble(WrappedEnsemble):
    """
    Alters trajectories given as arguments by taking Python slices.
    """

    def __init__(self, ensemble, region):
        super(SlicedTrajectoryEnsemble, self).__init__(ensemble)
        if type(region) == int:
            if region == -1:
                self.region = slice(region, None)
            else:
                self.region = slice(region, region + 1)
        else:
            self.region = region

    def _alter(self, trajectory):
        return trajectory[self.region]

    def _str(self):
        # TODO: someday may add different string support for slices with
        # only one frame
        start = "" if self.region.start is None else str(self.region.start)
        stop = "" if self.region.stop is None else str(self.region.stop)
        step = "" if self.region.step is None else " every " + str(
            self.region.step)
        return ("(" + str(self.ensemble) +
                " in {" + start + ":" + stop + "}" + step + ")")


class SuffixTrajectoryEnsemble(WrappedEnsemble):
    """
    Ensemble which prepends its trajectory to a given trajectory.

    Used in backward shooting.
    """

    def __init__(self, ensemble, add_trajectory):
        super(SuffixTrajectoryEnsemble, self).__init__(ensemble)
        self.add_trajectory = add_trajectory
        self._cached_trajectory = paths.Trajectory(add_trajectory.as_proxies())

    def _alter(self, trajectory):
        logger.debug("Starting Suffix._alter")
        # logger.debug(
        #     "altered " + str([id(i) for i in self._cached_trajectory]))
        reset = self._cache_can_prepend.check(trajectory)
        # logger.debug(
        #     "altered " + str([id(i) for i in self._cached_trajectory]))
        # logger.debug("traj    " + str([id(i) for i in trajectory]))
        # logger.debug("trajrev " + str([id(i) for i in trajectory.reversed]))
        # reset = False
        if not reset:
            logger.debug("SuffixTrajectory was not reset")
            first_frame = trajectory.get_as_proxy(-1)
            if self._cached_trajectory.get_as_proxy(0) != first_frame:
                self._cached_trajectory.insert(0, first_frame)
        else:
            self._cached_trajectory = trajectory.reversed + self.add_trajectory

        # logger.debug("revtraj " + str([id(i) for i in revtraj]))
        # logger.debug("add     " + str([id(i) for i in self.add_trajectory]))
        # logger.debug(
        #     "altered " + str([id(i) for i in self._cached_trajectory]))

        return self._cached_trajectory

    def can_append(self, trajectory, trusted=None):
        raise RuntimeError("SuffixTrajectoryEnsemble.can_append is nonsense.")

    def strict_can_append(self, trajectory, trusted=None):
        # was overridden in WrappedEnsemble: here should raise same error as
        # can_append does
        return self.can_append(trajectory, trusted)


class PrefixTrajectoryEnsemble(WrappedEnsemble):
    """
    Ensemble which appends its trajectory to a given trajectory.

    Used in forward shooting.
    """

    def __init__(self, ensemble, add_trajectory):
        super(PrefixTrajectoryEnsemble, self).__init__(ensemble)
        self.add_trajectory = add_trajectory
        self._cached_trajectory = paths.Trajectory(add_trajectory.as_proxies())

    def _alter(self, trajectory):
        logger.debug("Starting _alter")
        reset = self._cache_can_append.check(trajectory)
        if not reset:
            final_frame = trajectory.get_as_proxy(-1)
            if self._cached_trajectory.get_as_proxy(-1) != final_frame:
                self._cached_trajectory.append(final_frame)
        else:
            logger.debug("doing it oldstyle")
            self._cached_trajectory = self.add_trajectory + trajectory

            # DEBUG
            # logger.debug("add   " + str([i for i in self.add_trajectory]))
            # logger.debug("traj  " + str([i for i in trajectory]))
            # logger.debug("cache " + str([i for i in self._cached_trajectory]))
            # oldstyle = self.add_trajectory + trajectory
            # for (t,b) in zip(self._cached_trajectory,
            # self.add_trajectory+trajectory):
            # logger.debug(str(t) + " ?=? " + str(b))
            # assert(t == b)
        # assert(len(self._cached_trajectory) == len(oldstyle))

        return self._cached_trajectory

    def can_prepend(self, trajectory, trusted=None):
        raise RuntimeError("PrefixTrajectoryEnsemble.can_prepend is nonsense.")

    def strict_can_prepend(self, trajectory, trusted=None):
        # was overridden in WrappedEnsemble: here should raise same error as
        # can_append does
        return self.can_prepend(trajectory, trusted)


class ReversedTrajectoryEnsemble(WrappedEnsemble):
    """
    Ensemble based on reversing the trajectory.
    """

    def _alter(self, trajectory):
        return trajectory.reverse()


class AppendedNameEnsemble(WrappedEnsemble):
    """
    Add string to ensemble name: allows multiple copies of an ensemble.
    """

    def __init__(self, ensemble, label):
        self.label = label
        super(AppendedNameEnsemble, self).__init__(ensemble)

    def _str(self):
        return str(self.ensemble) + " " + self.label


class OptionalEnsemble(WrappedEnsemble):
    """
    An ensemble which is optional for SequentialEnsembles.
    """

    def __init__(self, ensemble):
        super(OptionalEnsemble, self).__init__(ensemble)
        self._new_ensemble = LengthEnsemble(0) | self.ensemble

    def _str(self):
        return "{" + str(self.ensemble) + "} (OPTIONAL)"


class SingleFrameEnsemble(WrappedEnsemble):
    """
    Convenience ensemble to `and` a LengthEnsemble(1) with a given ensemble.

    Frequently used for SequentialEnsembles.

    Attributes
    ----------
    ensemble : :class:`openpathsampling.ensemble.Ensemble`
        the ensemble which should be represented in the single frame

    Notes
    -----
    We allow the user to choose to be stupid: if, for example, the user
    tries to make a SingleFrameEnsemble from an ensemble which requires
    more than one frame to be satisfied (e.g., a SequentialEnsemble with
    more than one subensemble), it can be created, but no path will ever
    satisfy it. Since we can't stop all possible mistakes, we don't bother
    here.
    """

    def __init__(self, ensemble):
        super(SingleFrameEnsemble, self).__init__(ensemble)
        self._new_ensemble = LengthEnsemble(1) & self.ensemble

    def _str(self):
        return "{" + str(self.ensemble) + "} (SINGLE FRAME)"


class MinusInterfaceEnsemble(WrappedEnsemble):
    """
    This creates an ensemble for the minus interface.

    The specific implementation allows us to use the multiple-segment minus
    ensemble described by Swenson and Bolhuis. The minus interface was
    originally developed by van Erp. For more details, see the section
    "Anatomy of a PathMover: the Minus Move" in the OpenPathSampling
    Documentation.

    Parameters
    ----------
    state_vol : :class:`.Volume`
        The Volume which defines the state for this minus interface
    innermost_vols : list of :class:`.Volume`
        The Volume defining the innermost interface with which this minus
        interface does its replica exchange.
    n_l : integer (greater than one)
        The number of segments crossing innermost_vol for this interface.

    References
    ----------
    T.S. van Erp. Phys. Rev. Lett.

    D.W.H. Swenson and P.G. Bolhuis. J. Chem. Phys. 141, 044101 (2014).
    doi:10.1063/1.4890037
    """

    # don't store unnecessary stuff we recreate at initialization
    # TODO: Check with David if it makes sense to store these and allow
    # them being used in __init__ instead of the self-made ones

    def __init__(self, state_vol, innermost_vols, n_l=2, forbidden=None,
                 greedy=False):
        if n_l < 2:
            raise ValueError("The number of segments n_l must be at least 2")

        self.state_vol = state_vol
        try:
            innermost_vols = list(innermost_vols)
        except TypeError:
            innermost_vols = [innermost_vols]

        if forbidden is None:
            forbidden = [paths.EmptyVolume()]
        else:
            try:
                forbidden = list(forbidden)
            except TypeError:
                forbidden = [forbidden]

        self.forbidden = forbidden
        forbidden_volume = paths.join_volumes(forbidden)
        forbidden_ensemble = paths.AllOutXEnsemble(forbidden_volume)

        self.innermost_vols = innermost_vols
        self.innermost_vol = paths.FullVolume()
        for vol in self.innermost_vols:
            self.innermost_vol = self.innermost_vol & vol
        self.greedy = greedy
        in_A = AllInXEnsemble(state_vol)
        out_A = AllOutXEnsemble(state_vol)
        in_X = AllInXEnsemble(self.innermost_vol)
        leave_X = PartOutXEnsemble(self.innermost_vol)
        # interstitial = out_A & in_X
        interstitial = self.innermost_vol - state_vol
        in_interstitial = AllInXEnsemble(interstitial)
        segment_ensembles = [paths.TISEnsemble(state_vol, state_vol, inner)
                             for inner in self.innermost_vols]

        self._segment_ensemble = join_ensembles(segment_ensembles)

        # interstitial = AllInXEnsemble(self.innermost_vol - state_vol)
        start = [
            SingleFrameEnsemble(in_A),
            OptionalEnsemble(in_interstitial),
        ]
        loop = [
            out_A, # & leave_X, # redundant b/c next stop for previous
            in_X  # & hitA # redundant due to stop req for previous outA
        ]
        end = [
            out_A, #  & leave_X,
            OptionalEnsemble(in_interstitial),
            SingleFrameEnsemble(in_A)
        ]
        sequence = start + loop * (n_l - 1) + end

        ensemble = paths.SequentialEnsemble(sequence) & forbidden_ensemble

        self.n_l = n_l

        super(MinusInterfaceEnsemble, self).__init__(ensemble)

    def to_dict(self):
        dct = super(MinusInterfaceEnsemble, self).to_dict()
        dct['state_vol'] = self.state_vol
        dct['innermost_vols'] = self.innermost_vols
        dct['innermost_vol'] = self.innermost_vol
        dct['_segment_ensemble'] = self._segment_ensemble
        dct['forbidden'] = self.forbidden
        dct['n_l'] = self.n_l
        return dct

    @property
    def extendable_sub_ensembles(self):
        # A-X-A and the one from TISEnsemble

        state_vol = self.state_vol

        sub_ensembles = {}

        in_A = AllInXEnsemble(state_vol)
        out_A = AllOutXEnsemble(state_vol)

        # this code is for potential

        # in_X = AllInXEnsemble(self.innermost_vol)
        # leave_X = PartOutXEnsemble(self.innermost_vol)
        # interstitial = out_A & in_X
        # segment_ensembles = [paths.TISEnsemble(state_vol, state_vol, inner)
        #                      for inner in self.innermost_vols]

        # start = [
        #     SingleFrameEnsemble(in_A),
        #     OptionalEnsemble(interstitial),
        # ]
        # loop = [
        #     out_A & leave_X,
        #     in_X  # & hitA # redundant due to stop req for previous outA
        # ]
        # end = [
        #     out_A & leave_X,
        #     OptionalEnsemble(interstitial),
        #     SingleFrameEnsemble(in_A)
        # ]

        # do not add higher orders, you would
        # for n_l in range(self.n_l - 2, 0, -1):
        #     # add ens with less loops
        #     sub_ensembles.append(
        #         SequentialEnsemble(start + loop * n_l + end))

        sub_ensembles['complex'] = self._segment_ensemble

        # and the simplest possible just crossing from in_state to outside
        sub_ensembles['minimal'] = \
            LengthEnsemble(2) &  \
            SequentialEnsemble([
                SingleFrameEnsemble(in_A),
                SingleFrameEnsemble(out_A)
            ])

        return sub_ensembles

    # def populate_minus_ensemble(self, partial_traj, minus_replica_id, engine):
    #     """
    #     Generate a sample for the minus ensemble by extending `partial_traj`
    #
    #     Parameters
    #     ----------
    #     partial_traj : :class:`openpathsampling.trajectory.Trajectory`
    #         trajectory to extend
    #     minus_replica_id : int or str
    #         replica ID for this sample
    #     engine : :class:`openpathsampling.dynamicsengine.DynamicsEngine`
    #         engine to use for MD extension
    #     """
    #     last_frame = partial_traj[-1]
    #     if not self._segment_ensemble(partial_traj):
    #         raise RuntimeError(
    #             "Invalid input trajectory for minus extension. (Not A-to-A?)"
    #         )
    #     fwd_extend_ens = PrefixTrajectoryEnsemble(self, partial_traj)
    #     extension = engine.generate(last_frame,
    #                                 [fwd_extend_ens.can_append])
    #     first_minus = paths.Trajectory(partial_traj + extension[1:])
    #     assert self(first_minus)
    #     minus_samp = paths.Sample(
    #         replica=minus_replica_id,
    #         trajectory=first_minus,
    #         ensemble=self
    #     )
    #     logger.info(first_minus.summarize_by_volumes_str(
    #         {"A": self.state_vol,
    #          "I": ~self.state_vol & self.innermost_vol,
    #          "X": ~self.innermost_vol})
    #     )
    #     return minus_samp

    # def populate_minus_ensemble_from_set(self, samples, minus_replica_id,
    #                                      engine):
    #     """
    #     Generate a sample for this minus ensemble by extending trajectory.
    #
    #     Parameters
    #     ----------
    #     samples : iterable of :class:`.Sample`
    #         samples with trajectories that might be extended
    #     minus_replica_id : int or str
    #         replica ID for the return sample
    #     engine : :class:`openpathsampling.dynamicsengine.DynamicsEngine`
    #         engine to use for MD extension
    #
    #     Returns
    #     -------
    #     :class:`.Sample` :
    #         a sample for this minus ensemble
    #     """
    #     partials = [s.trajectory for s in samples
    #                 if self._segment_ensemble(s.trajectory)]
    #     if len(partials) == 0:
    #         # TODO: add support for trying to run backwards
    #         raise RuntimeError("No trajectories can be extended")
    #
    #     samp = None
    #
    #     good_sample = False
    #     while not good_sample:
    #         partial_traj = partials[0]
    #         # I think it should be impossible to RuntimeError in this
    #         samp = self.populate_minus_ensemble(
    #             partial_traj=partial_traj,
    #             minus_replica_id=minus_replica_id,
    #             engine=engine
    #         )
    #
    #         good_sample = samp.ensemble(samp.trajectory)
    #
    #     return samp


class TISEnsemble(WrappedEnsemble):
    """An ensemble for TIS (or AMS).

    Begin in `initial_states`, end in either `initial_states` or
    `final_states`, and cross `interface`.

    Attributes
    ----------
    initial_states : `openpathsampling.volume.Volume` or list of `openpathsampling.volume.Volume`
        Volume(s) that only the first or last frame may be in
    final_states : `openpathsampling.volume.Volume` or list of `openpathsampling.volume.Volume`
        Volume(s) that only the last frame may be in
    interface : `openpathsampling.volume.Volume`
        Volume which the trajectory must exit to be accepted
    orderparameter : `openpathsampling.collectivevariable.CollectiveVariable`
        CV to be used as order parameter for this
    """

    @property
    def extendable_sub_ensembles(self):
        # this is tricky. The only extendable sub-ensembles are (In, Out)
        # at the crossing of leaving or entering the core

        # pick only the initial ones like for A-X-AB pick A
        states = list(set(self.initial_states))

        volume = paths.volume.join_volumes(states)

        return {
            'minimal':
            LengthEnsemble(2) &
            SequentialEnsemble([
                SingleFrameEnsemble(AllInXEnsemble(volume)),
                SingleFrameEnsemble(AllOutXEnsemble(volume))
            ])
        }

    def __init__(self, initial_states, final_states, interface,
                 orderparameter=None, cv_max=None, lambda_i=None):
        # regularize to list of volumes
        # without orderparameter, some info can't be obtained
        try:
            _ = len(initial_states)
        except TypeError:
            initial_states = [initial_states]

        try:
            _ = len(final_states)
        except TypeError:
            final_states = [final_states]

        volume_a = paths.volume.join_volumes(initial_states)
        volume_b = paths.volume.join_volumes(final_states)

        ensemble = SequentialEnsemble([
            AllInXEnsemble(volume_a) & LengthEnsemble(1),
            OptionalEnsemble(AllOutXEnsemble(volume_a | volume_b)),
            AllInXEnsemble(volume_a | volume_b) & LengthEnsemble(1)
        ]) & PartOutXEnsemble(interface)
        super(TISEnsemble, self).__init__(ensemble)

        self.initial_states = initial_states
        self.final_states = final_states
        self.interface = interface
        #        self.name = interface.name
        self.orderparameter = orderparameter  # TODO: is this used? remove?
        self.cv_max = cv_max
        self.lambda_i = lambda_i
        self._initial_volumes = volume_a
        self._final_volumes = volume_b | volume_a

    def __call__(self, trajectory, trusted=None, candidate=False):
        logger.debug("TIS ENSEMBLE: candidate={0}".format(str(candidate)))
        use_candidate = (candidate and self.lambda_i is not None)
        if use_candidate and self.cv_max is not None:
            logger.debug("Using candidate shortcut with self.cv_max")
            return (
                self._initial_volumes(trajectory[0])
                & self._final_volumes(trajectory[-1])
                & (self.cv_max(trajectory) > self.lambda_i)
            )
        elif use_candidate and self.orderparameter is not None:
            logger.debug("Using candidate shortcut with max(orderparameter)")
            # as a candidate trajectory, we assume that only the first and
            # final frames can be in a state
            #logger.debug("initial: " +
                         #str(self._initial_volumes(trajectory[0])))
            #logger.debug("final: " +
                         #str(self._final_volumes(trajectory[0])))
            #logger.debug("max: " +
                         #str(max(self.orderparameter(trajectory))))
            return (
                self._initial_volumes(trajectory[0])
                & self._final_volumes(trajectory[-1])
                & (max(self.orderparameter(trajectory)) > self.lambda_i)
            )
        else:
            logger.debug("No shortcut possible")
            # it still works fine if we use the slower algorithm
            return super(TISEnsemble, self).__call__(trajectory, trusted)

    def trajectory_summary(self, trajectory):
        initial_state_i = None
        final_state_i = None
        for state_i in range(len(self.initial_states)):
            if self.initial_states[state_i](trajectory.get_as_proxy(0)):
                initial_state_i = state_i
                break
        all_states = self.initial_states + self.final_states
        for state_i in range(len(all_states)):
            if all_states[state_i](trajectory.get_as_proxy(-1)):
                final_state_i = state_i
                break

        if self.orderparameter is not None:
            lambda_traj = self.orderparameter(trajectory)
            min_lambda = min(lambda_traj)
            max_lambda = max(lambda_traj)
        else:
            min_lambda = None
            max_lambda = None

        return {
            'initial_state': initial_state_i,
            'final_state': final_state_i,
            'max_lambda': max_lambda,
            'min_lambda': min_lambda
        }

    def trajectory_summary_str(self, trajectory):
        summ = self.trajectory_summary(trajectory)
        all_states = self.initial_states + self.final_states
        # TODO: remove the .name from this when string returns correctly
        init_st_i = summ['initial_state']
        fin_st_i = summ['final_state']
        # TODO: how can we have None?
        if init_st_i is None:
            init_st = "None"
        else:
            init_st = str(self.initial_states[summ['initial_state']].name)
        if fin_st_i is None:
            fin_st = "None"
        else:
            fin_st = str(all_states[summ['final_state']].name)

        # if self.orderparameter is not None:
        #     opname = self.orderparameter.name
        # else:
        #     opname = "None"
        min_l = str(summ['min_lambda'])
        max_l = str(summ['max_lambda'])
        mystr = (
            "initial_state=" + init_st + " " +
            "final_state=" + fin_st + " " +
            "min_lambda=" + min_l + " " +
            "max_lambda=" + max_l + " "
        )
        return mystr

    def _str(self):
        return str(self.ensemble)


# class EnsembleFactory(object):
    # """
    # Convenience class to construct Ensembles
    # """

    # @staticmethod
    # def StartXEnsemble(volume):
        # """
        # Construct an ensemble that starts (x[0]) in the specified volume

        # Parameters
        # ----------
        # volume : :class:`openpathsampling.volume.Volume`
            # The volume to start in

        # Returns
        # -------
        # ensemble : :class:`openpathsampling.ensemble.Ensemble`
            # The constructed Ensemble
        # """
        # return AllInXEnsemble(volume, 0)

    # @staticmethod
    # def EndXEnsemble(volume):
        # """
        # Construct an ensemble that ends (x[-1]) in the specified volume

        # Parameters
        # ----------
        # volume : :class:`openpathsampling.volume.Volume`
            # The volume to end in

        # Returns
        # -------
        # ensemble : :class:`openpathsampling.ensemble.Ensemble`
            # The constructed Ensemble
        # """
        # return AllInXEnsemble(volume, -1)

    # @staticmethod
    # def A2BEnsemble(volume_a, volume_b, trusted=True):
        # """
        # Construct an ensemble that starts in `volume_a`, ends in
        # `volume_b` and is in either volumes in between

        # Parameters
        # ----------
        # volume_a : :class:`openpathsampling.Volume`
            # The volume to start in
        # volume_b : :class:`openpathsampling.Volume`
            # The volume to end in

        # Returns
        # -------
        # ensemble : :class:`openpathsampling.Ensemble`
            # The constructed Ensemble
        # """
        # # TODO: this is actually only for flexible path length TPS now
        # return SequentialEnsemble([
            # SingleFrameEnsemble(AllInXEnsemble(volume_a)),
            # AllOutXEnsemble(volume_a | volume_b),
            # SingleFrameEnsemble(AllInXEnsemble(volume_b))
        # ])

    # @staticmethod
    # def TISEnsembleSet(volume_a, volume_b, volumes_x, orderparameter,
                       # lambdas=None):
        # if lambdas is None:
            # lambdas = [None] * len(volumes_x)
        # myset = [paths.TISEnsemble(volume_a, volume_b, vol, orderparameter,
                                   # lambda_i)
                 # for (vol, lambda_i) in zip(volumes_x, lambdas)]
        # return myset
