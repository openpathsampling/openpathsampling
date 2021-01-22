from openpathsampling.netcdfplus import StorableNamedObject
from openpathsampling.engines.dynamics_engine import DynamicsEngine
from openpathsampling.engines.snapshot import BaseSnapshot, SnapshotDescriptor
from openpathsampling.engines.toy import ToySnapshot
import numpy as np
import os

import logging

import psutil
import signal
import shlex
import time

import sys
if sys.version_info > (3, ):
    long = int

import linecache

logger = logging.getLogger(__name__)

def close_file_descriptors(basename):
    """Close file descriptors for the given filename.

    There may be instances when the file reading leaves a file open (like if
    an error was encountered while reading the file). This closes any open
    file descriptors.
    """
    len_name = len(basename)
    fds = [p for p in psutil.Process().open_files()
           if p.path[-len_name:] == basename]
    for fd in fds:
        logger.debug("Closing " + fd.path)
        open(fd.fd).close()


def _debug_open_files(where=None, ext=".trr"):
    len_ext = len(ext)
    if logger.isEnabledFor(logging.DEBUG):
        logger.debug("REACHED " + where)
        open_trrs = [p.path for p in psutil.Process().open_files()
                     if ext == p.path[-len_ext:]]
        if len(open_trrs) > 0:
            message = []
            for open_trr in open_trrs:
                loc_msg = "Open file in {}: {}".format(str(where), open_trr)
                message.append(loc_msg)

            raise Exception("\n".join(message))

def _debug_snapshot_loading(snapshot):
    snapshot.load_details()
    snapshot.clear_cache()


class FilenameSetter(StorableNamedObject):
    """Just use numbers, as we did previously.

    This is the default for compatibility reasons, but it not recommended.
    Generally, we recommend using :class:`.RandomString`.
    """
    # the weird use of this as the base class is because engine options has
    # some weird type testing that requires replace object to be instances
    # of the default
    def __init__(self, count=0):
        super(FilenameSetter, self).__init__()
        self.count = count

    def __call__(self):
        retval = self.count
        self.count += 1
        return retval

    def reset(self, count=0):
        self.count = count


class RandomStringFilenames(FilenameSetter):
    """Use a random string for filename prefixes.

    This is recommended, and will become the default in future versions of
    OpenPathSampling.

    Parameters
    ----------
    length : int
        number of character in the resulting string
    allowed : str
        string containing the allowed characters to return
    """
    _allowed = 'abcdefghijklmnopqrstuvwxyz0123456789'
    def __init__(self, length=8, allowed=None):
        super(RandomStringFilenames, self).__init__()
        self.length = length
        allowed = allowed if allowed is not None else self._allowed
        self.allowed = np.array([a for a in allowed])

    def __call__(self):
        return "".join(np.random.choice(self.allowed, self.length))

class _InternalizedEngineProxy(DynamicsEngine):
    """Wrapper that allows snapshots to be "internalized."

    This is needed because the dynamically-registered snapshot tables are
    associated to an engine: each engine creates exactly one snapshot table.
    We use the internalized engine to create those snapshots.
    """
    def __init__(self, engine):
        descriptor = SnapshotDescriptor.construct(
            snapshot_class=engine.InternalizedSnapshotClass,
            snapshot_dimensions=engine.descriptor._dimensions
        )
        self.engine = engine
        super(_InternalizedEngineProxy, self).__init__(options={},
                                                       descriptor=descriptor)

    def to_dict(self):
        # We only need the engine to load again
        dct = {'engine': self.engine}
        return dct

    @property
    def name(self):
        return self.engine.name + " (internalized)"

    def __getattr__(self, attr):
        return getattr(self.engine, attr)


class ExternalEngine(DynamicsEngine):
    """
    Generic object to handle arbitrary external engines. Subclass to use.
    """

    _default_options = {
        'n_frames_max' : 10000,
        'name_prefix' : "test",
        'default_sleep_ms' : 100,
        'auto_optimize_sleep' : True,
        'engine_sleep' : 100,
        'engine_directory' : "",
        'n_spatial' : 1,
        'n_atoms' : 1,
        'n_poll_per_step': 1,
        'filename_setter': FilenameSetter(),
    }

    killsig = signal.SIGTERM

    def __init__(self, options, descriptor, template,
                 first_frame_in_file=False):
        # needs to be overridden for each engine
        super(ExternalEngine, self).__init__(options=options,
                                             descriptor=descriptor)
        self.template = template
        self.sleep_ms = self.default_sleep_ms
        self.start_time = None
        self.first_frame_in_file = first_frame_in_file
        self._traj_num = -1
        self._current_snapshot = template
        self.n_frames_since_start = None
        self.internalized_engine = _InternalizedEngineProxy(self)

    def to_dict(self):
        dct = super(ExternalEngine, self).to_dict()
        # this is a trick to avoid circular references: the idea is that we
        # create an identical internalized engine in our __init__ (i.e., the
        # engine's from_dict) -- therefore we create it from scratch each
        # time and just save/set the UUID
        dct.update({
            'internalized_engine_uuid': str(self.internalized_engine.__uuid__)
        })
        return dct

    @classmethod
    def from_dict(cls, dct):
        dct = dct.copy()
        internalized_engine_uuid = dct.pop('internalized_engine_uuid', None)
        obj = super(ExternalEngine, cls).from_dict(dct)
        if internalized_engine_uuid:
            # TODO: in future versions of OPS, get_/set_uuid should always
            # return strings
            obj.internalized_engine.__uuid__ = long(internalized_engine_uuid)
        return obj

    @property
    def current_snapshot(self):
        return self._current_snapshot

    @current_snapshot.setter
    def current_snapshot(self, snap):
        self._current_snapshot = snap

    def generate_next_frame(self):
        # should be completely general
        next_frame_found = False
        logger.debug("Looking for frame %d", self.n_frames_since_start+1)
        while not next_frame_found:
            try:
                next_frame = self.read_frame_from_file(self.output_file,
                                                       self.frame_num)
            except IOError:
                # maybe the file doesn't exist
                if self.proc.is_running():
                    logger.info("Waiting for file to be written")
                    next_frame = None
                else:
                    raise
            #print self.frame_num, next_frame # DEBUG LOGGER
            now = time.time()
            if next_frame == "partial":
                if self.proc.poll() is not None:
                    raise RuntimeError("External engine died unexpectedly")
                time.sleep(0.001) # wait a millisec and rerun
            elif next_frame is None:
                if self.proc.poll() is not None:
                    raise RuntimeError("External engine died unexpectedly")
                logger.debug("Sleeping for {:.2f}ms".format(self.sleep_ms))
                time.sleep(self.sleep_ms/1000.0)
            elif isinstance(next_frame, BaseSnapshot): # success
                self.n_frames_since_start += 1
                logger.debug("Found frame %d", self.n_frames_since_start)
                self.current_snapshot = next_frame
                next_frame_found = True
                self.frame_num += 1
            else:  # pragma: no cover
                raise RuntimeError("Strange return value from read_next_frame_from_file")
            if self.auto_optimize_sleep and self.n_frames_since_start > 0:
                n_poll_per_step = self.options['n_poll_per_step']
                elapsed = now - self.start_time
                time_per_step = elapsed / self.n_frames_since_start
                self.sleep_ms = time_per_step / n_poll_per_step * 1000.0

        return self.current_snapshot

    def start(self, snapshot=None):
        super(ExternalEngine, self).start(snapshot)
        self._traj_num += 1
        self.frame_num = 0
        self.n_frames_since_start = 0
        file_prefix = self.filename_setter()
        self.set_filenames(file_prefix)
        self.write_frame_to_file(self.input_file, self.current_snapshot, "w")
        self.prepare()

        cmd = shlex.split(self.engine_command())
        self.start_time = time.time()
        try:
            logger.info(self.engine_command())
            # TODO: add the ability to have handlers for stdin and stdout
            self.proc = psutil.Popen(shlex.split(self.engine_command()),
                                     preexec_fn=os.setsid)
        except OSError:  # pragma: no cover
            raise  #TODO: need to handle this, but do what?
        else:
            logger.info("Started engine: " + str(self.proc))

        if self.first_frame_in_file:
            _ = self.generate_next_frame()  # throw away repeat first frame

    def stop(self, trajectory):
        super(ExternalEngine, self).stop(trajectory)
        logger.info("total_time {:.4f}".format(time.time() - self.start_time))
        proc = self.who_to_kill()
        logger.info("About to send signal %s to %s", str(self.killsig),
                    str(proc))
        try:
            proc.send_signal(self.killsig)
            logger.debug("Signal has been sent")
            proc.wait()  # wait for the zombie to die
            logger.debug("Zombie should be dead")
        except psutil.NoSuchProcess:
            logger.debug("Tried to kill process, but it was already dead")
        self.cleanup()

    # FROM HERE ARE THE FUNCTIONS TO OVERRIDE IN SUBCLASSES:
    def read_frame_from_file(self, filename, frame_num):
        """Reads given frame number from file, and returns snapshot.

        If no frame is available, returns None. If the frame appears to be
        partially written, returns string "partial".
        """
        raise NotImplementedError()

    def write_frame_to_file(self, filename, snapshot, mode="a"):
        """Writes given snapshot to file."""
        raise NotImplementedError()

    def who_to_kill(self):
        """Returns psutil.Process object to send kill signal to.

        Might override to send kill signal to a process other than the one
        directly spawned above (e.g., when launching parallel runs)
        """
        # this should only be called if you're about to kill the process; if
        # the process doesn't exist, you shouldn't be killing anything and
        # it will raise an error
        return self.proc

    def prepare(self):
        """
        Any preparation between writing snapshot and running command
        """
        pass

    def cleanup(self):
        """Any cleanup actions to do after the subprocess dies."""
        pass

    def set_filenames(self, number):
        """Sets names for files associated with trajectory `number`"""
        # pass instead of not implemented because maybe you don't need it
        # for some engine? usually you do, though
        pass

    def engine_command(self):
        """Generates a string for the command to run the engine."""
        raise NotImplementedError()

