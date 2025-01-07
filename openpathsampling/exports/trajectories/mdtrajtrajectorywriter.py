from .core import TrajectoryWriter
from openpathsampling.integration_tools import HAS_MDTRAJ
from collections.abc import Iterable


class MDTrajTrajectoryWriter(TrajectoryWriter):
    """Generic save to MDTraj.

    Note that this will not include velocities, and therefore isn't suitable
    for saving data that could be used in a restart.
    """
    def __init__(self, mdtraj_selection=None):
        if not HAS_MDTRAJ:  # -no-cov-
            raise ImportError("MDTraj is not available")

        self.mdtraj_selection = mdtraj_selection
        self._sel = None  # the actual atom indices of selection

    def _write(self, trajectory, filename):
        mdt = trajectory.to_mdtraj()

        if self._sel is None:
            if self.mdtraj_selection is None:
                pass
            elif isinstance(self.mdtraj_selection, str):
                self._sel = mdt.topology.select(self.mdtraj_selection)
            elif isinstance(self.mdtraj_selection, Iterable):
                self._sel = list(self.mdtraj_selection)
            else:
                raise TypeError("mdtraj_selection must be a string or an "
                                f"iterable; got '{self.mdtraj_selection}'")

        if self._sel is not None:
            mdt = mdt.atom_slice(self._sel)

        mdt.save(filename)
