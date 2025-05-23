"""
Gromacs support in OpenPathSampling

OpenPathSampling supports Gromacs as an "external engine," meaning that
Gromacs itself runs as a coprocess, and OPS uses the file system as an
intermediary to monitor when a path needs to be terminated.
"""

import logging
logger = logging.getLogger(__name__)

from openpathsampling.integration_tools import md, HAS_MDTRAJ
if HAS_MDTRAJ:
    TRRTrajectoryFile = md.formats.TRRTrajectoryFile

from openpathsampling.engines import ExternalEngine
from openpathsampling.engines import features
from openpathsampling.engines.snapshot import BaseSnapshot, SnapshotDescriptor
from openpathsampling.engines.topology import MDTrajTopology
from openpathsampling.engines.external_snapshots import \
        ExternalMDSnapshot, InternalizedMDSnapshot
from openpathsampling.tools import ensure_file
from openpathsampling.exports.trajectories import TRRTrajectoryWriter

import os
import psutil
import shlex
import time
import numpy as np

from openpathsampling.engines.external_engine import (
    _debug_open_files, close_file_descriptors
)

def _remove_file_if_exists(filename):  # pragma: no cover
    #  not requiring coverage here because it's part of Gromacs integration;
    #  gets covered if gmx is installed though
    if os.path.isfile(filename):
        os.remove(filename)

class _GroFileEngine(ExternalEngine):
    SnapshotClass = ExternalMDSnapshot
    InternalizedSnapshotClass = InternalizedMDSnapshot
    def __init__(self, gro):
        self.gro = gro
        traj = md.load(gro)
        self.topology = MDTrajTopology(traj.topology)
        n_atoms = self.topology.n_atoms
        n_spatial = self.topology.n_spatial
        descriptor = SnapshotDescriptor.construct(
            snapshot_class=self.SnapshotClass,
            snapshot_dimensions={'n_spatial': n_spatial,
                                 'n_atoms': n_atoms}
        )
        super(_GroFileEngine, self).__init__(options={},
                                            descriptor=descriptor,
                                            template=None)

    def to_dict(self):
        return {'gro': self.gro}

    @classmethod
    def from_dict(cls, dct):
        return cls(dct['gro'])

    def read_frame_data(self, file_name, file_position):
        traj = md.load(file_name)
        xyz = traj.xyz[0]
        vel = np.zeros(shape=xyz.shape)
        box = traj.unitcell_vectors[0]
        return (xyz, vel, box)



def snapshot_from_gro(gro_file):
    template_engine = _GroFileEngine(gro_file)
    snapshot = ExternalMDSnapshot(file_name=gro_file,
                                  file_position=0,
                                  engine=template_engine)
    return snapshot


class GromacsEngine(ExternalEngine):
    """
    External engine wrapper for Gromacs (using indirect API).

    This provides Gromacs support, using our indirect engine API (TODO
    link).

    OPS runs Gromacs based on the mdp file that you provide. This mdp file
    MUST output in the TRR format, and the velocity and position save
    frequency in the TRR must be the same (that is, you need to have
    ``nstxout`` = ``nstvout``, and they must not be 0).

    Parameters
    ----------
    gro : string
        File for the grompp ``-c`` flag. This is often a .gro, but note that
        you may get better support with other integrations (e.g., MDTraj) if
        you use a PDB.
    mdp : string
        .mdp file
    top : string
        .top file
    options : dict
        Dictionary of option name to value. Gromacs-specific option names
        are

            * ``gmx_executable``: Prefix to gromacs commands, which are run
              as ``{gmx_prefix}command``. Default is 'gmx ' (note the
              space).  This allows you to use either Gromacs 4 or Gromacs 5,
              as well as specifying the path to your version of Gromacs.
            * ``grompp_args``: Additional arguments to ``grompp``. The
              defaults take ``-c {self.gro} -f {self.mdp} -p {self.top}
              -t {self.input_file}``, where the input filename is set by
              :meth:`.set_filenames`. Default is the empty string.
            * ``mdrun_args``: Additional arguments to ``mdrun``. The
              defaults take ``-s topol.top -o self.output_file
              -e self.edr_file -g self.log_file``, where the ``topol.top``
              is generated by :meth:`.prepare`, and the other filenames are
              set by :meth:`.set_filenames`. Default is the empty string.
            * ``snapshot_timestep``: time between frames analysed by ops.
              You keep track of the unit, I'd advise ps so the output
              rates will be in ps. Example. 2 fs timestep in the mdp with
              nstxout of 30 would give snapshot_timestep of 60 fs = 0.06 ps

    base_dir : string
        root directory where all files will be found (defaults to pwd)
    prefix : string
        prefix within ``base_dir`` for output folders (defaults to gmx)
    """
    _default_options = dict(ExternalEngine._default_options,
        **{
            'gmx_executable': "gmx ",
            'grompp_args': "",
            'mdrun_args': "",
            'snapshot_timestep':1.0
        }
    )
    GROMPP_CMD = ("{e.options[gmx_executable]}grompp -c {e.gro} "
                  + "-f {e.mdp} -p {e.top} -t {e.input_file} "
                  + "-po {e.mdout_file} -o {e.tpr_file} "
                  + "{e.options[grompp_args]}")
    MDRUN_CMD = ("{e.options[gmx_executable]}mdrun -s {e.tpr_file} "
                 + "-o {e.output_file} -e {e.edr_file} -g {e.log_file} "
                 + "{mdrun_args}")
    # use these as CMD.format(e=engine, **engine.options)
    SnapshotClass = ExternalMDSnapshot
    InternalizedSnapshotClass = InternalizedMDSnapshot
    clear_snapshot_cache = True
    def __init__(self, gro, mdp, top, options, base_dir="", prefix="gmx"):
        self.base_dir = base_dir
        self.gro = os.path.join(base_dir, gro)
        self.mdp = os.path.join(base_dir, mdp)
        self.top = os.path.join(base_dir, top)
        self.prefix = os.path.join(base_dir, prefix)

        self.gro_contents, self._gro_hash = ensure_file(self.gro)
        self.mdp_contents, self._mdp_hash = ensure_file(self.mdp)
        self.top_contents, self._top_hash = ensure_file(self.top)

        # TODO: move to a later stage, before first traj
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

        # initial placeholders
        self.input_file = "INITIAL.trr"
        self.output_file = self.prefix + "_trr/OUTPUT_NAME.trr"
        self.edr_file = self.prefix + "_edr/OUTPUT_NAME.edr"
        self.log_file = self.prefix + "_log/OUTPUT_NAME.log"
        self.tpr_file = "topol.tpr"
        self.mdout_file = "mdout.mdp"

        self._mdtraj_topology = None

        super(GromacsEngine, self).__init__(options, descriptor, template,
                                             first_frame_in_file=True)

    def _default_trajectory_writer(self):
        return TRRTrajectoryWriter()

    def to_dict(self):
        dct = super(GromacsEngine, self).to_dict()
        local_dct = {
            'gro': self.gro,
            'mdp': self.mdp,
            'top': self.top,
            'options': self.options,
            'base_dir': self.base_dir,
            'prefix': self.prefix,
            'gro_contents': self.gro_contents,
            'mdp_contents': self.mdp_contents,
            'top_contents': self.top_contents,
        }
        dct.update(local_dct)
        return dct

    @classmethod
    def from_dict(cls, dct):
        dct = dict(dct)  # make a copy
        for ftype in ['gro', 'top', 'mdp']:
            contents = dct.pop(ftype + "_contents")
            _ = ensure_file(filename=dct[ftype], old_contents=contents)
        return super(GromacsEngine, cls).from_dict(dct)

    @property
    def mdtraj_topology(self):
        if self._mdtraj_topology:
            return self._mdtraj_topology
        return self.topology.mdtraj

    @mdtraj_topology.setter
    def mdtraj_topology(self, value):
        self._mdtraj_topology = value

    def read_frame_data(self, filename, frame_num):
        """
        Returns pos, vel, box or raises error
        """
        # if self._last_filename != filename:
            # try:
                # self._file.close()
            # except AttributeError:
                # pass  # first time thru, self._file is None
            # self._file = TRRTrajectoryFile(filename)
        # f = self._file
        # do we need to reopen the TRR each time to avoid problems with the
        # fiel length changing?
        trr = TRRTrajectoryFile(filename)
        f = trr
        logger.debug("Reading file %s frame %d (of %d)",
                     filename, frame_num, len(f))
        # logger.debug("File length: %d", len(f))
        try:
            f.seek(offset=frame_num)
            data = f._read(n_frames=1, atom_indices=None, get_velocities=True)
        finally:
            trr.close()

        return data[0][0], data[5][0], data[3][0]

    def read_frame_from_file(self, file_name, frame_num):
        # note: this only needs to return the file pointers -- but should
        # only do so once that frame has been written!
        basename = os.path.basename(file_name)
        # basename should be in the format [0-9]+\.trr (as set by the
        # trajectory_filename method)
        # file_number = int(basename.split('.')[0])
        try:
            xyz, vel, box = self.read_frame_data(file_name, frame_num)
        except (IndexError, OSError, IOError) as e:
            # this means that no such frame exists yet, so we return None
            # IndexError in older version, OSError more recently (specific
            # MDTraj error)
            logger.debug("Expected exception caught: " + str(e))
            close_file_descriptors(basename)
            return None
        except RuntimeError as e:
            # TODO: matches "TRR read error"
            logger.debug("Received partial frame for %s %d", file_name,
                         frame_num+1)
            return 'partial'
        else:
            logger.debug("Creating snapshot")
            snapshot =  ExternalMDSnapshot(file_name=file_name,
                                           file_position=frame_num,
                                           engine=self)
            return snapshot

    def write_frame_to_file(self, filename, snapshot, mode='w'):
        if os.path.isfile(filename):
            # stop if we already have this file; could also happen because
            # of a weird behavior in a mover. You must remove the files if
            # you don't want them.
            raise RuntimeError("File " + str(filename) + " exists. "
                               + "Preventing overwrite.")
        # type control before passing things to Cython code
        xyz = np.asarray([snapshot.xyz], dtype=np.float32)
        time = np.asarray([0.0], dtype=np.float32)
        step = np.asarray([0], dtype=np.int32)
        box = np.asarray([np.asarray(snapshot.box_vectors)], dtype=np.float32)
        lambd = np.asarray([0.0], dtype=np.float32)
        vel = np.asarray([np.asarray(snapshot.velocities)], dtype=np.float32)

        try:
            trr = TRRTrajectoryFile(filename, mode)
            trr._write(xyz, time, step, box, lambd, vel)
        finally:
            trr.close()
            close_file_descriptors(filename)

    def trajectory_filename(self, number):
        # TODO: remove code path allowing ints (API break for 2.0)
        trr_dir = self.prefix + "_trr/"
        if isinstance(number, int):
            num_str = num_str = '{:07d}'.format(number)
        else:
            num_str = number
        return trr_dir + num_str + '.trr'

    def set_filenames(self, number):
        if isinstance(number, int):
            num_str = '{:07d}'.format(number + 1)
            self.output_file = self.trajectory_filename(number + 1)
            init_filename = "initial_frame.trr"
            self.filename_setter.reset(number + 1)
        else:
            num_str = number
            self.output_file = self.trajectory_filename(num_str)
            init_filename = num_str + "_initial_frame.trr"
            self.mdout = num_str + "_mdout.mdp"
            self.tpr_file = num_str + "_" + "topol.tpr"

        self.input_file = os.path.join(self.base_dir, init_filename)
        self.edr_file = os.path.join(self.prefix + "_edr", num_str + '.edr')
        self.log_file = os.path.join(self.prefix + "_log", num_str + '.log')

    @property
    def grompp_command(self):
        cmd = self.GROMPP_CMD.format(e=self)
        return cmd

    def prepare(self):  # pragma: no cover
        # coverage ignored b/c Travis won't have gmx. However, we do have a
        # test that covers this if gmx is present (otherwise it is skipped)

        # ensure that the files haven't changed before we run trajectory
        _ = ensure_file(self.gro, self.gro_contents, self._gro_hash)
        _ = ensure_file(self.mdp, self.mdp_contents, self._mdp_hash)
        _ = ensure_file(self.top, self.top_contents, self._top_hash)

        # grompp and mdrun
        cmd = self.grompp_command
        logger.info(cmd)
        run_cmd = shlex.split(cmd)
        return_code = psutil.Popen(run_cmd, preexec_fn=os.setsid).wait()
        return return_code

    def cleanup(self):  # pragma: no cover
        # tested when traj is run, which we don't on CI
        _remove_file_if_exists(self.input_file)
        _remove_file_if_exists(self.tpr_file)
        _remove_file_if_exists(self.mdout_file)

    def engine_command(self):
        args = self.options['mdrun_args'].format(prev_traj=self._traj_num-1,
                                                 next_traj=self._traj_num)
        cmd = self.MDRUN_CMD.format(e=self, mdrun_args=args)
        return cmd
