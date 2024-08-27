from collections import namedtuple
import logging
import openpathsampling as paths

from .options import create_default_options, canonicalize_mover
from .mpl import ThinLines
from .svg import SVGLines

Connector = namedtuple("Connector", ['frame', 'mccycle'])
PathTreeStep = namedtuple(
    "PathTreeStep",
    ['offset', 'mccycle', 'mover', 'trajectory', 'accepted',
     'connector', 'segments']
)
PathTreeStep.__doc__ = """Step information needed for path tree.

Each :class:`.PathTreeStep` becomes a row of the path tree.

Parameters
----------
offset : int
    the left edge of where the trajectory will be drawn
mccycle : int
    the MC step number for this step
mover : :class:`.PathMover`
    the path mover associated with this step
trajectory : :class:`.Trajectory`
    the trajectory that was generated in this step
accepted : bool
    whether the step was accepted
connector : :class:`.Connector`
    information connecting the the shooting betwee this step and its
    original step
segments: :class:`.TrajectorySegments`

"""
TrajectorySegments = namedtuple("TrajectorySegments", ['existing', 'new'])


def make_layout_blocks(path_tree_steps):
    """Separate these into blocks separated by None offsets.

    If no frames are in common between successive steps, we create a new
    "block". Each block will be centered independently.
    """
    steps = iter(path_tree_steps)
    first_step = next(steps)
    minimum = 0
    maximum = len(first_step.trajectory)
    blocks = []
    block = [first_step]
    for step in steps:
        if step.offset is None:
            block.append(block, maximum - minimum)
            minumum = 0
            maximum = len(step.trajectory)
            block = []
        else:
            minimum = min(minimum, offset)
            maximum = max(maxumum, offset + len(step.trajectory))

        block.append(step)

    return blocks

def make_segments(prev, trial):
    """Identify new and existing segments.
    """
    prev_traj = set(prev)  # as set for lookup performance
    start = None
    prev_mode = None
    results = {'existing': [], 'new': []}
    for idx, snap in enumerate(trial):
        mode = {True: 'existing', False: 'new'}[snap in prev_traj]
        if mode != prev_mode:
            if start is not None:
                results[prev_mode].append((start, idx))
            start = idx
        prev_mode = mode

    if prev_mode is not None:
        results[prev_mode].append((start, idx + 1))

    return TrajectorySegments(**results)

def _shooting_points(details):
    try:
        shooting_pt = details.shooting_snapshot
    except AttributeError:
        # not a shooting move
        return None, None

    try:
        # two-way shooting
        mod_shooting = details.modified_shooting_snapshot
    except AttributeError:
        # one-way shooting
        mod_shooting = shooting_pt

    return shooting_pt, mod_shooting

def _offset_from_shooting(prev, trial, shooting_pt, mod_shooting):
    prev_idx = prev.trajectory.index(shooting_pt)
    new_idx = trial.index(mod_shooting)
    return  prev_idx - new_idx

def _offset_from_trajectories(prev, trial):
    shared = prev.trajectory.shared_configurations(trial,
                                                   time_reversal=True)
    if not shared:
        return None

    # I'm not entirely sure that this approach will always work
    min_prev = min(prev.trajectory.index(s) for s in shared)
    min_trial = min(trial.index(s) for s in shared)
    offset = min_prev - min_trial
    return offset

def _replica_or_ensemble(replica, ensemble):
        if replica is None and ensemble is None:
            replica = 0
        if replica is not None and ensemble is not None:
            raise ValueError("Can't give both replica and ensemble")
        if replica is None:
            selector = ensemble
        else:
            selector = replica
        return selector


class PathTree(paths.progress.SimpleProgress):
    def __init__(self, steps, replica=None, ensemble=None):
        self.selector = _replica_or_ensemble(replica, ensemble)
        if steps is not None:
            self.path_tree_steps = self.generate_path_tree_steps(steps)

    def _filter_steps(self, steps):
        for step in steps:
            trials = paths.SampleSet(step.change.trials)
            try:
                trial = trials[self.selector].trajectory
            except KeyError:
                continue
            yield trial, step

    @staticmethod
    def _get_shift(prev, trial, shooting_pt, mod_shooting):
        if shooting_pt is None:
            # not a shooting move
            offset = _offset_from_trajectories(prev, trial)
        else:
            # shooting move
            offset = _offset_from_shooting(prev, trial, shooting_pt,
                                           mod_shooting)

        return offset

    def _first_path_tree_step(self, step):
        """Generate the tree_step for the first step"""
        trajectory = step.active[self.selector].trajectory
        return PathTreeStep(
            offset=0,
            mccycle=step.mccycle,
            mover=None,
            accepted=True,
            connector=None,
            trajectory=trajectory,
            segments=TrajectorySegments(existing=[],
                                        new=[(0, len(trajectory))]),
        )

    def _connector(self, shooting_pt, mod_shooting, trial, existing):
        if shooting_pt is None:
            return None

        connector_idx = trial.index(mod_shooting)
        mccycle = existing[shooting_pt]
        return Connector(connector_idx, mccycle)

    def generate_path_tree_steps(self, steps):
        steps = iter(self.progress(steps))
        prev = self._first_path_tree_step(next(steps))
        path_tree_steps = [prev]
        connectors = []
        pre_existing = {snap: prev.mccycle for snap in prev.trajectory}
        offset = 0
        for trial, step in self._filter_steps(steps):
            # track when the snapshot first entered (continuously)
            current = {snap: step.mccycle for snap in trial}
            existing = dict(current)
            existing.update(pre_existing)

            # get information to make the tree step
            details = step.change.canonical.details
            shooting_pt, mod_shooting = _shooting_points(details)

            shift = self._get_shift(prev, trial, shooting_pt, mod_shooting)
            connector = self._connector(shooting_pt, mod_shooting, trial,
                                        existing)

            if shift is None:
                offset = 0
                shift = 0

            tree_step = PathTreeStep(
                offset=offset + shift,
                mccycle=step.mccycle,
                mover=step.change.canonical.mover,
                accepted=step.change.accepted,
                trajectory=trial,
                connector=connector,
                segments=make_segments(prev.trajectory, trial)
            )

            if step.change.accepted:
                offset += shift

            path_tree_steps.append(tree_step)

            if tree_step.accepted:
                prev = tree_step
                pre_existing = {k: existing[k] for k in current}

        return path_tree_steps


    def svg(self, plotter=None):
        if plotter is None:
            plotter = SVGLines

        if type(plotter) is type:
            plotter = plotter()

        return plotter.draw(self)

    def matplotlib(self, plotter=None):
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise RuntimeError("matplotlib is not installed")

        if plotter is None:
            plotter = ThinLines

        if type(plotter) is type:
            plotter = plotter()

        ax = plotter.ax
        plotter.draw_trajectories(self.path_tree_steps)

        return ax.figure
