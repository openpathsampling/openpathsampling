"""
Gromacs support in OpenPathSampling

OpenPathSampling supports Gromacs as an "external engine," meaning that
Gromacs itself runs as a coprocess, and OPS uses the file system as an
intermediary to monitor when a path needs to be terminated.
"""

import logging

from mdtraj.formats import TRRTrajectoryFile

from openpathsampling.engines import ExternalEngine
from openpathsampling.engines import features
from openpathsampling.engines.snapshot import BaseSnapshot
import features as gmx_features

# TODO: all gmx_features should be moved to external_md

@features.base.attach_features([
    features.engine,
    gmx_features.coordinates,
    # gmx_features.velocities,
    # gmx_features.box_vectors
    gmx_features.file_info
])
class GromacsSnapshot(BaseSnapshot):
    """
    Snapshot for external MD engines
    """
    def __init__(self, file_number=None, file_position=None):
        self.file_number = file_number
        self.file_position = file_position
        self._xyz = None
        self._velocities = None
        self._box_vectors = None

    def load_details(self)
        filename = self.engine.trajectory_filename(self.file_number)
        (xyz, vel, box) = self.engine.read_frame_data(filename,
                                                      self.file_position)
        self._xyz = xyz
        self._velocities = vel
        self._box_vectors = box

    def clear_cache(self):
        self._xyz = None
        self._velocities = None
        self._box_vectors = None


class Gromacs5Engine(ExternalEngine):
    def __init__(self, gro, mdp, options, template, name="gmx"):
        self.gro = gro
        self.mdp = mdp
        self.name = name
        self._last_file = None  # file open/close efficiency
        self._last_file_name = None
        # TODO: update options with correct n_spatial, n_atoms
        # TODO: add snapshot_timestep; first via options, later read mdp
        super(Gromacs5Engine, self).__init__(options, template)

    def read_frame_data(self, filename, frame_num):
        """
        Returns pos, vel, box or raises error
        """
        if self._last_filename != filename:
            self._last_file.close()
            self._file = TRRTrajectoryFile(file_name)
        self._file.seek(offset=frame_num)
        data = self._file._read(n_frames=1, atom_indices=None,
                                get_velocities=True)
        return data[0], data[5], data[3]

    def read_frame_from_file(self, filename, frame_num):
        # note: this only needs to return the file pointers -- but should
        # only do so once that frame has been written!
        file_number = filename  # TODO: get actual file name
        try:
            self.read_frame_data(filename, frame_num)
        except Exception as e:
            # TODO: what kind of exception is this?
            print e
            return 'partial'
        else:
            return GromacsSnapshot(file_number, file_position)

    def write_frame_to_file(self, filename, snapshot, mode='a'):
        pass

    def trajectory_filename(self, number):
        trr_dir = self.name + "_trr/"
        return trr_dir + '{:07d}'.format(number) + '.trr'

    def set_filenames(self, number):
        pass

    def engine_command(self):
        pass
