"""
@author: JD Chodera
@author: JH Prinz
"""

import numpy as np

from openpathsampling.integration_tools import (
    error_if_no_mdtraj, is_simtk_quantity_type, md
)
from openpathsampling.netcdfplus import StorableObject, LoaderProxy
import openpathsampling as paths


# ==============================================================================
# TRAJECTORY
# ==============================================================================


class Trajectory(list, StorableObject):
    """
    Simulation trajectory. Essentially a python list of snapshots
    """

    engine = None

    def __init__(self, trajectory=None):
        """
        Create a simulation trajectory object

        Parameters
        ----------

        trajectory : :obj:`Trajectory` or list of :obj:`openpathsampling.engines.BaseSnapshot`
            if specified, make a deep copy of specified trajectory
        """

        # Initialize list.
        list.__init__(self)
        StorableObject.__init__(self)

        if trajectory is not None:
            if type(trajectory) is Trajectory:
                self.extend(trajectory.iter_proxies())
            else:
                self.extend(trajectory)

    def extend(self, iterable):
        if type(iterable) is Trajectory:
            list.extend(self, iterable.iter_proxies())
        else:
            list.extend(self, iterable)

    def to_dict(self):
        return {
            'snapshots': self.as_proxies()
        }

    @classmethod
    def from_dict(cls, dct):
        return cls(dct['snapshots'])

    def __str__(self):
        return 'Trajectory[' + str(len(self)) + ']'

    def __repr__(self):
        return 'Trajectory[' + str(len(self)) + ']'

    @property
    def snapshots(self):
        return list(self)

    def map(self, fnc, allow_fast=True):
        """
        This runs a function and tries to be fast.

        Fast here means that functions that are purely based on CVs can be
        evaluated without actually loading the real Snapshot object. This
        functions tries to do that and if it fails it does it the usual way
        and creates the snapshot object. This bears the possibility that
        the function uses the fake snapshots and returns a non-sense value.
        It is up to the user to make sure this will not happen.
        """

        if allow_fast:
            try:
                return [fnc(frame) for frame in list.__iter__(self)]
            except:
                pass

        return [fnc(frame) for frame in self]

    @property
    def reversed(self):
        """
        Returns a reversed (shallow) copy of the trajectory itself. Effectively
        creates a new Trajectory object and then fills it with shallow reversed
        copies of the contained snapshots.

        Returns
        -------
        :class:`openpathsampling.trajectory.Trajectory`
            the reversed trajectory
        """

        return Trajectory([snap for snap in reversed(self)])

    @property
    def n_snapshots(self):
        """
        Return the number of frames in the trajectory.

        Returns
        -------
        length (int) - the number of frames in the trajectory

        Notes
        -----
        Might be removed in later versions for len(trajectory) is more pythonic

        See also
        --------
        len

        """

        return len(self)

    def __getattr__(self, item):
        """
        Fallback to access Snapshot properties

        """
        if len(self) == 0:
            return []

        snapshot_class = self[0].__class__
        def is_snapshot_attr(cls, item):
            return hasattr(cls, item) or (hasattr(cls, '__features__') and
                                         item in cls.__features__.variables)

        if is_snapshot_attr(snapshot_class, item):
            # if there's a trajectory_item in features, that should be a
            # function to return a trajectory
            traj_item = "trajectory_" + item
            if traj_item in snapshot_class.__features__.functions:
                traj_func = getattr(snapshot_class, traj_item)
                return traj_func(self)

            # get the results
            out = [getattr(snap, item) for snap in self]

            # if the first result is a numpy object, return the whole as a
            # numpy array
            if isinstance(out[0], np.ndarray):
                out = np.array(out)

            return out

        # behavior when it can't be delegated to snapshot
        std_msg = "'{0}' object has no attribute '{1}'"
        snap_msg = "Cannot delegate to snapshots. "
        snap_msg += "'{2}' has no attribute '{1}'"
        spacer = "\n                "
        msg = (std_msg + spacer + snap_msg).format(
            str(self.__class__.__name__),
            item,
            snapshot_class.__name__
        )
        raise AttributeError(msg)


    # ==========================================================================
    # LIST INHERITANCE FUNCTIONS
    # ==========================================================================

    def __getslice__(self, *args, **kwargs):
        ret = list.__getslice__(self, *args, **kwargs)
        if type(ret) is list:
            ret = Trajectory(ret)

        return ret

    # this is intuitive. hash(Trajectory(traj)) == hash(traj)
    # but hash(LoaderProxy(..., traj.__uuid__)) != hash(traj)

    def __hash__(self):
        if len(self) == 0:
            return hash(tuple())
        else:
            return hash(
                (list.__getitem__(self, 0), len(self),
                 list.__getitem__(self, -1)))

    # this might be faster, but does not allow to compare arbitrary
    # trajectories as one might expect. hash(Trajectory(traj)) != hash(traj)
    # but hash(LoaderProxy(..., traj.__uuid__)) == hash(traj)
    # it could also lead to better caching and memory behaviour

    # __hash__ = StorableObject.__hash__
    #
    # __eq__ = StorableObject.__eq__
    # __ne__ = StorableObject.__ne__
    #
    def __getitem__(self, index):
        # Allow for numpy style selection using lists
        if hasattr(index, '__iter__'):
            ret = [list.__getitem__(self, i) for i in index]
        else:
            ret = list.__getitem__(self, index)

        if type(ret) is list:
            ret = Trajectory(ret)
        elif type(ret) is LoaderProxy:
            ret = ret.__subject__

        return ret

    def __reversed__(self):
        for snap_idx in range(len(self) - 1, -1, -1):
            yield self[snap_idx].reversed

    def index_symmetric(self, value):
        """
        Return index of a snapshot or its reversed inside a trajectory

        """
        try:
            fw = self.index(value)
        except ValueError:
            fw = None
            pass

        try:
            bw = self.index(value.reversed)
        except ValueError:
            bw = None
            pass

        if fw is None:
            if bw is None:
                raise KeyError(
                    '%r or its reversed is not found in trajectory.')
            else:
                return bw
        else:
            if bw is None:
                return fw
            else:
                return min(fw, bw)

    def contains_symmetric(self, item):
        """
        Test whether a snapshot or its reversed is in a trajectory

        Returns
        -------
        bool

        """
        fw = item in self
        if not fw:
            return item.reversed in self
        else:
            return True

    def get_as_proxy(self, item):
        """
        Get an actual contained element

        This will also return lazy proxy objects and not the referenced ones
        as does __iter__, __reversed__ or __getitem__. Useful for faster access
        to the elements

        This is equal to use list.__getitem__(trajectory, item)

        Returns
        -------
        :obj:`Snapshot` or :obj:`openpathsampling.netcdfplus.proxy.LoaderProxy`
        """
        return list.__getitem__(self, item)

    def as_proxies(self):
        """
        Returns all contains all actual elements

        This will also return lazy proxy objects and not the references ones
        as does __iter__, __reversed__ or __getitme__. Useful for faster access
        to the elements

        Returns
        -------
        list of :obj:`Snapshot` or :obj:`openpathsampling.netcdfplus.LoaderProxy`

        """
        return list(self.iter_proxies())

    def iter_proxies(self):
        """
        Returns an iterator over all actual elements

        This will also return lazy proxy objects and not the references ones
        as does __iter__, __reversed__ or __getitme__. Useful for faster
        access to the elements

        Returns
        -------
        Iterator() over list of :class:`openpathsampling.snapshot.Snapshot`
        or :class:`openpathsampling.netcdfplus.proxy.LoaderProxy`


        """
        return list.__iter__(self)

    def __iter__(self):
        """
        Return an iterator over all snapshots in the storage

        This will always give real :class:`openpathsampling.snapshot.Snapshot`
        objects and never proxies to snapshots.
        If you prefer proxies (if available) use `.items()`

        Returns
        -------
        Iterator()
            The iterator that iterates the objects in the store

        """

        for snap_idx in range(len(self)):
            yield self[snap_idx]

    def __add__(self, other):
        t = Trajectory(self)
        t.extend(other)
        return t

    # ==========================================================================
    # PATH ENSEMBLE FUNCTIONS
    # ==========================================================================

    def summarize_by_volumes(self, label_dict):
        """Summarize trajectory based on number of continuous frames in volumes.

        This uses a dictionary of disjoint volumes: the volumes must be disjoint
        so that every frame can be mapped to one volume. If the frame maps to
        none of the given volumes, it returns the label None.

        Parameters
        ----------
        label_dict : dict
            dictionary with labels for keys and volumes for values

        Returns
        -------
        list of tuple
            format is (label, number_of_frames)
        """
        last_vol = None
        count = 0
        segment_labels = []
        # list.__iter__ for speed
        for frame in list.__iter__(self):
            in_state = []
            for key in label_dict.keys():
                vol = label_dict[key]
                if vol(frame):
                    in_state.append(key)
            if len(in_state) > 1:
                raise RuntimeError(
                    "Volumes given to summarize_by_volumes not disjoint")
            if len(in_state) == 0:
                current_vol = None
            else:
                current_vol = in_state[0]

            if last_vol == current_vol:
                count += 1
            else:
                if count > 0:
                    segment_labels.append( (last_vol, count) )
                last_vol = current_vol
                count = 1
        segment_labels.append( (last_vol, count) )
        return segment_labels

    def summarize_by_volumes_str(self, label_dict, delimiter="-"):
        """
        Return string version of the volumes visited by this trajectory.

        See `Trajectory.summarize_by_volumes` for details.

        Parameters
        ----------
        label_dict : dict
            dictionary with labels for keys and volumes for values
        delimiter : string (default "-")
            string used to separate volumes in output

        Returns
        -------
        string
            order in which this trajectory visits the volumes in
            `label_dict`, separated by the `delimiter`
        """
        summary = self.summarize_by_volumes(label_dict)
        return delimiter.join([str(s[0]) for s in summary])

    # ==========================================================================
    # ANALYSIS FUNCTIONS
    # ==========================================================================

    def is_correlated(self, other, time_reversal=False):
        """
        Checks if two trajectories share a common snapshot

        Parameters
        ----------
        other : :class:`openpathsampling.trajectory.Trajectory`
            the second trajectory to check for common snapshots

        Returns
        -------
        bool
            returns True if at least one snapshot appears in both trajectories
        """
        return bool(self.shared_configurations(
            other,
            time_reversal=time_reversal))

    def shared_configurations(self, other, time_reversal=False):
        """
        Returns a set of shared snapshots

        Parameters
        ----------
        other : :class:`openpathsampling.trajectory.Trajectory`
            the second trajectory to use

        Returns
        -------
        set of :class:`openpathsampling.snapshot.Snapshot`
            the set of common snapshots
        """
        if not time_reversal:
            return set(self.as_proxies()) & set(other.as_proxies())
        else:
            return set(self.as_proxies()) & \
                (set(other.as_proxies()) |
                 set([snap.reversed for snap in other.as_proxies()]))

    def shared_subtrajectory(self, other, time_reversal=False):
        """
        Returns a subtrajectory which only contains frames present in other

        Parameters
        ----------
        other : :class:`openpathsampling.trajectory.Trajectory`
            the second trajectory to use

        Returns
        -------
        :class:`openpathsampling.trajectory.Trajectory`
            the shared subtrajectory
        """
        shared = self.shared_configurations(other, time_reversal=time_reversal)
        return Trajectory([snap for snap in self.iter_proxies() if snap in shared])

    def unique_subtrajectory(self, other):
        """
        Returns a subtrajectory which contains frames not present in other


        Parameters
        ----------
        other : :class:`openpathsampling.trajectory.Trajectory`
            the second trajectory to use

        Returns
        -------
        :class:`openpathsampling.trajectory.Trajectory`
            the unique frames subtrajectory (opposite of shared)
        """
        unique = set([snap for snap in self]) - set([snap for snap in other])
        return Trajectory([snap for snap in self if snap in unique])


    def subtrajectory_indices(self, subtrajectories):
        """
        Returns a list of lists of indices for frames from subtrajectories.

        Parameters
        ----------
        subtrajectories : list of :class:`.Trajectory`
            input list of subtrajectories

        Returns
        -------
        list of list of int
            the indices within this trajectory of the frames in each
            subtrajectory
        """
        if isinstance(subtrajectories, Trajectory):
            return [self.index(s) for s in subtrajectories]
        else:
            return [[self.index(s) for s in subtrj]
                    for subtrj in subtrajectories]

    # ==========================================================================
    # UTILITY FUNCTIONS
    # ==========================================================================

    def to_mdtraj(self, topology=None):
        """
        Construct a mdtraj.Trajectory object from the Trajectory itself

        Parameters
        ----------
        topology : :class:`mdtraj.Topology`
            If not None this topology will be used to construct the mdtraj
            objects otherwise the topology object will be taken from the
            configurations in the trajectory snapshots.

        Returns
        -------
        :class:`mdtraj.Trajectory`
            the trajectory

        Notes
        -----

        If the OPS trajectory is zero-length (has no snapshots), then this
        fails. OPS cannot currently convert zero-length trajectories to
        MDTraj, because an OPS zero-length trajectory cannot determine its
        MDTraj topology.
        """
        error_if_no_mdtraj("Converting to mdtraj")
        try:
            snap = self[0]
        except IndexError:
            raise ValueError("Cannot convert zero-length trajectory "
                             + "to MDTraj")
        if topology is None:
            # TODO: maybe add better error output?
            # if AttributeError here, engine doesn't support mdtraj
            topology = snap.engine.mdtraj_topology

        output = self.xyz

        traj = md.Trajectory(output, topology)
        box_vectors = self.box_vectors
        # box_vectors is a list with an entry for each frame of the traj
        # if they're all None, we return None, not [None, None, ..., None]
        if not np.any(box_vectors):
            box_vectors = None

        traj.unitcell_vectors = box_vectors
        return traj

    @property
    def topology(self):
        """
        Return a Topology object representing the topology of the
        current view of the trajectory

        Returns
        -------
        :class:`openpathsampling.topology.Topology`
            the topology object

        """

        topology = None

        if len(self) > 0 and self[0].topology is not None:
            # if no topology is defined
            topology = self[0].topology

        return topology

    @staticmethod
    def _to_list_of_trajectories(trajectories):
        if isinstance(trajectories, Trajectory):
            trajectories = [trajectories]
        elif isinstance(trajectories, paths.Sample):
            trajectories = [trajectories.trajectory]
        elif isinstance(trajectories, paths.SampleSet):
            trajectories = [s.trajectory for s in trajectories]
        elif isinstance(trajectories, list):
            if len(trajectories) > 0:
                trajectories = [
                    obj.trajectory if isinstance(obj, paths.Sample) else obj
                    for obj in trajectories
                    ]
        elif isinstance(trajectories, paths.BaseSnapshot):
            return paths.Trajectory([trajectories])
        elif isinstance(trajectories, paths.BaseSnapshot):
            return paths.Trajectory([trajectories])

        return trajectories
