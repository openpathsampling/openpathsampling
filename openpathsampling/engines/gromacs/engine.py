"""
Gromacs support in OpenPathSampling

OpenPathSampling supports Gromacs as an "external engine," meaning that
Gromacs itself runs as a coprocess, and OPS uses the file system as an
intermediary to monitor when a path needs to be terminated.
"""

import logging
logger = logging.getLogger(__name__)

from mdtraj.formats import TRRTrajectoryFile
import mdtraj as md

from openpathsampling.engines import ExternalEngine
from openpathsampling.engines import features
from openpathsampling.engines.snapshot import BaseSnapshot, SnapshotDescriptor
from openpathsampling.engines.openmm.topology import MDTrajTopology
from . import features as gmx_features

import os
import psutil
import shlex
import numpy as np

# TODO: all gmx_features should be moved to external_md; along with the
# snapshot

@features.base.attach_features([
    features.engine,
    gmx_features.coordinates,
    gmx_features.velocities,
    gmx_features.box_vectors,
    gmx_features.file_info
])
class ExternalMDSnapshot(BaseSnapshot):
    """
    Snapshot for external MD engines

    Internally, this only stores the file_name and the file_position. All
    specific details (positions, velocities, box vectors) are loaded from
    file when requested.

    Parameters
    ----------
    file_name : str
        the name of the external file where the positions/velocities/etc.
        reside
    file_position : int
        position within the file; the engine should be able to load data for
        this specific snapshot based on this number
    engine : :class:`.DynamicsEngine`
        the engine associated with this snapshot
    """
    def __init__(self, file_name=None, file_position=None, engine=None):
        # these are done in place of calling super
        self._reversed = None
        self.__uuid__ = self.get_uuid()
        # these are the requried attributes
        self.file_name = file_name
        self.file_position = file_position
        self.engine = engine
        self.velocity_direction = 1  # by default; reversed flips it
        # these are containers for temporary data
        self._xyz = None
        self._velocities = None
        self._box_vectors = None

    def load_details(self):
        """Cache coords, velocities, box vectors from the external file"""
        (xyz, vel, box) = self.engine.read_frame_data(self.file_name,
                                                      self.file_position)
        self._xyz = xyz
        self._velocities = vel
        self._box_vectors = box

    def set_details(self, xyz, velocities, box_vectors):
        """Set coords, velocities, and box vectors.

        This is mainly used if OPS must modify/create a snapshot.

        Parameters
        ----------
        xyz : np.array
            unitless coordinates
        velocities : np.array
            velocities
        box_vectors : np.array
            unit cell for the periodic box
        """
        try:
            self.load_details()
        except:
            pass
        else:
            raise RuntimeError("Can't set details if frame already exists.")
        finally:
            self._xyz = xyz
            self._velocities = velocities
            self._box_vectors = box_vectors

    def clear_cache(self):
        """Remove internal details from snapshot.

        These details should always be accessible later using
        :method:`.load_details`. Removing them allows them memory to be
        freed.
        """
        self._xyz = None
        self._velocities = None
        self._box_vectors = None

    def __repr__(self):
        num_str = "file_name=" + str(self.file_name)
        pos_str = "file_position=" + str(self.file_position)
        eng_str = "engine=" + repr(self.engine)
        args = ", ".join([num_str, pos_str, eng_str])
        return "{cls_str}(".format(cls_str=self.cls) + args + ")"


def snapshot_from_gro(gro_file):
    class GroFileEngine(ExternalEngine):
        def __init__(self, gro):
            traj = md.load(gro)
            self.topology = MDTrajTopology(traj.topology)
            n_atoms = self.topology.n_atoms
            n_spatial = self.topology.n_spatial
            descriptor = SnapshotDescriptor.construct(
                snapshot_class=ExternalMDSnapshot,
                snapshot_dimensions={'n_spatial': n_spatial,
                                     'n_atoms': n_atoms}
            )
            super(GroFileEngine, self).__init__(options={},
                                                descriptor=descriptor,
                                                template=None)

        def read_frame_data(self, file_name, file_position):
            traj = md.load(file_name)
            xyz = traj.xyz[0]
            vel = np.zeros(shape=xyz.shape)
            box = traj.unitcell_vectors[0]
            return (xyz, vel, box)

    template_engine = GroFileEngine(gro_file)
    snapshot = ExternalMDSnapshot(file_name=gro_file,
                                  file_position=0,
                                  engine=template_engine)
    return snapshot


class GromacsEngine(ExternalEngine):
    """
    External engine wrapper for Gromacs (using indirect API).

    This provides Gromacs support, using our indirect engine API (TODO
    link).

    Parameters
    ----------
    gro : string
        .gro file
    mdp : string
        .mdp file
    top : string
        .top file
    base_dir : string
        root directory where all files will be found (defaults to pwd)
    options : dict
        Dictionary of option name to value. Gromacs-specific option names
        are
            * ``gmx_prefix``: Prefix to gromacs commands, which are run as
              ``{gmx_prefix}command``. Default is 'gmx ' (note the space).
              This allows you to use either Gromacs 4 or Gromacs 5, as well
              as specifying the path to your version of Gromacs.
            * ``grompp_args``: Additional arguments to ``grompp``. The
              defaults take ``-c {self.gro} -f {self.mdp} -p {self.top}
              -t {self.input_file}``, where the input filename is set by
              :meth:`.set_filenames`. Default is the empty string.
            * ``mdrun_args``: Additional arguments to ``mdrun``. The
              defaults take ``-s topol.top -o self.output_file
              -e self.edr_file -g self.log_file``, where the ``topol.top``
              is generated by :meth:`.prepare`, and the other filenames are
              set by :meth:`.set_filenames`. Default is the empty string.
    """
    _default_options = dict(ExternalEngine._default_options,
        **{
            'gmx_executable': "gmx ",
            'grompp_args': "",
            'mdrun_args': ""
        }
    )
    def __init__(self, gro, mdp, top, options, base_dir="", prefix="gmx"):
        self.base_dir = base_dir
        self.gro = os.path.join(base_dir, gro)
        self.mdp = os.path.join(base_dir, mdp)
        self.top = os.path.join(base_dir, top)
        self.prefix = os.path.join(base_dir, prefix)

        dirs = [self.prefix + s for s in ['_trr', '_log', '_edr']]
        for d in dirs:
            try:
                os.mkdir(d)
            except OSError:
                pass  # the directory already exists

        self._file = None  # file open/close efficiency
        self._last_filename = None
        # TODO: add snapshot_timestep; first via options, later read mdp
        template = snapshot_from_gro(self.gro)
        self.topology = template.topology
        descriptor = template.engine.descriptor  # descriptor from gro file

        super(GromacsEngine, self).__init__(options, descriptor, template,
                                             first_frame_in_file=True)

    def to_dict(self):
        # TODO: eventually, read files in and store the text
        return {
            'base_dir': self.base_dir,
            'gro': self.gro,
            'mdp': self.mdp,
            'top': self.top,
            'options': self.options,
            'base_dir': self.base_dir,
            'prefix': self.prefix
        }

    def read_frame_data(self, filename, frame_num):
        """
        Returns pos, vel, box or raises error
        """
        if self._last_filename != filename:
            try:
                self._file.close()
            except AttributeError:
                pass  # first time thru, self._file is None
            self._file = TRRTrajectoryFile(filename)
        self._file.seek(offset=frame_num)
        data = self._file._read(n_frames=1, atom_indices=None,
                                get_velocities=True)
        return data[0][0], data[5][0], data[3][0]

    def read_frame_from_file(self, file_name, frame_num):
        # note: this only needs to return the file pointers -- but should
        # only do so once that frame has been written!
        basename = os.path.basename(file_name)
        # basename should be in the format [0-9]+\.trr (as set by the
        # trajectory_filename method)
        file_number = int(basename.split('.')[0])
        try:
            xyz, vel, box = self.read_frame_data(file_name, frame_num)
        except (IndexError, OSError):
            # this means that no such frame exists yet, so we return None
            # IndexError in older version, OSError more recently (specific
            # MDTraj error)
            return None
        except RuntimeError as e:
            # TODO: matches "TRR read error"
            return 'partial'
        else:
            return ExternalMDSnapshot(file_name=file_name,
                                      file_position=frame_num,
                                      engine=self)

    def write_frame_to_file(self, filename, snapshot, mode='w'):
        if os.path.isfile(filename):
            # stop if we already have this file; could also happen because
            # of a weird behavior in a mover. You must remove the files if
            # you don't want them.
            raise RuntimeError("File " + str(filename) + " exists. "
                               + "Preventing overwrite.")
        trr = TRRTrajectoryFile(filename, mode)
        # type control before passing things to Cython code
        xyz = np.asarray([snapshot.xyz], dtype=np.float32)
        time = np.asarray([0.0], dtype=np.float32)
        step = np.asarray([0], dtype=np.int32)
        box = np.asarray([snapshot.box_vectors], dtype=np.float32)
        lambd = np.asarray([0.0], dtype=np.float32)
        vel = np.asarray([snapshot.velocities], dtype=np.float32)
        trr._write(xyz, time, step, box, lambd, vel)
        trr.close()

    def trajectory_filename(self, number):
        trr_dir = self.prefix + "_trr/"
        return trr_dir + '{:07d}'.format(number) + '.trr'

    def set_filenames(self, number):
        self.input_file = os.path.join(self.base_dir, "initial_frame.trr")
        self.output_file = self.trajectory_filename(number + 1)
        num_str = '{:07d}'.format(number + 1)
        self.edr_file = os.path.join(self.prefix + "_edr", num_str + '.edr')
        self.log_file = os.path.join(self.prefix + "_log", num_str + '.log')

    @property
    def grompp_command(self):
        cmd = "{gmx}grompp -c {gro} -f {mdp} -p {top} -t {inp} {xtra}".format(
            gmx=self.options['gmx_executable'], gro=self.gro, mdp=self.mdp,
            top=self.top, inp=self.input_file,
            xtra=self.options['grompp_args']
        )
        return cmd

    def prepare(self):  # pragma: no cover
        # coverage ignored b/c Travis won't have gmx. However, we do have a
        # test that covers this if gmx is present (otherwise it is skipped)
        cmd = self.grompp_command
        logger.info(cmd)
        run_cmd = shlex.split(cmd)
        return_code = psutil.Popen(run_cmd, preexec_fn=os.setsid).wait()
        return return_code

    def cleanup(self):
        if os.path.isfile(self.input_file):
            os.remove(self.input_file)

    def engine_command(self):
        # gmx mdrun -s topol.tpr -o trr/0000001.trr -g 0000001.log
        args = self.options['mdrun_args'].format(prev_traj=self._traj_num-1,
                                                 next_traj=self._traj_num)
        cmd = "{gmx}mdrun -s topol.tpr -o {out} -e {edr} -g {log} {args}"
        cmd = cmd.format(gmx=self.options['gmx_executable'],
                         out=self.output_file,
                         edr=self.edr_file,
                         log=self.log_file,
                         args=args)
        return cmd
