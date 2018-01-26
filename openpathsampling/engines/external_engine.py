from openpathsampling.engines.dynamics_engine import DynamicsEngine
from openpathsampling.engines.snapshot import BaseSnapshot
from openpathsampling.engines.toy import ToySnapshot
import numpy as np
import os

import logging

import psutil
import signal
import shlex
import time

import linecache

import sys # DEBUG

logger = logging.getLogger(__name__)

class ExternalEngine(DynamicsEngine):
    """
    Generic object to handle arbitrary external engines.

    Typically, this will be subclassed for any given engine. As written, it
    will work with the trivial `engine.c` developed for testing purposes.
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
        'n_poll_per_step': 1
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

    @property
    def current_snapshot(self):
        return self._current_snapshot

    @current_snapshot.setter
    def current_snapshot(self, snap):
        self._current_snapshot = snap

    def generate_next_frame(self):
        # should be completely general
        next_frame_found = False
        logger.debug("Looking for frame")
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
                time.sleep(0.001) # wait a millisec and rerun
            elif next_frame is None:
                if not self.proc.is_running():
                    raise RuntimeError("External engine died unexpectedly")
                logger.info("Sleeping for {:.2f}ms".format(self.sleep_ms))
                time.sleep(self.sleep_ms/1000.0)
            elif isinstance(next_frame, BaseSnapshot): # success
                self.n_frames_since_start += 1
                logger.info("Found frame")
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
        self.set_filenames(self._traj_num)
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

        if self.first_frame_in_file:
            _ = self.generate_next_frame()  # throw away repeat first frame

    def stop(self, trajectory):
        super(ExternalEngine, self).stop(trajectory)
        logger.info("total_time {:.4f}".format(time.time() - self.start_time))
        proc = self.who_to_kill()
        proc.send_signal(self.killsig)
        proc.wait()  # wait for the zombie to die
        self.cleanup()

    # FROM HERE ARE THE FUNCTIONS TO OVERRIDE IN SUBCLASSES:
    def read_frame_from_file(self, filename, frame_num):
        """Reads given frame number from file, and returns snapshot.

        If no frame is available, returns None. If the frame appears to be
        partially written, returns string "partial".
        """
        # under most circumstances, start with linecache.checkcache and
        # setting the value of the first line
        linecache.checkcache(filename)
        first_line = frame_num + 1

        # create a snapshot out of lines starting with first_line... if
        # nothing exists, linecache returns '', so we return None.
        # Otherwise, try to make a snapshot and return "partial" if we fail
        line = linecache.getline(filename, first_line)
        if line is '':
            snap = None
        else:
            try:
                splitted = line.split()
                if len(splitted) == 2:
                    coords = float(splitted[0])
                    vels = float(splitted[1])
                else:
                    raise ValueError()  # force the raise we then ignore
                snap = ToySnapshot(coordinates=np.array([[coords]]),
                                   velocities=np.array([[vels]]))
            except ValueError:
                snap = "partial"
        return snap

    def write_frame_to_file(self, filename, snapshot, mode="a"):
        """Writes given snapshot to file."""
        f = open(filename, mode)
        snapshot_text = "{pos} {vel}\n".format(pos=snapshot.xyz[0][0],
                                               vel=snapshot.velocities[0][0])
        f.write(snapshot_text)
        f.close()

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
        self.input_file = self.name_prefix + str(number) + ".inp"
        self.output_file = self.name_prefix + str(number) + ".out"

    def engine_command(self):
        """Generates a string for the command to run the engine."""
        if self.engine_directory != "":
            engine_path = os.path.join(self.engine_directory, "engine")
        else:  # pragma: no cover
            engine_path = "engine"
        return (engine_path + " " + str(self.engine_sleep)
                + " " + str(self.output_file) + " " + str(self.input_file))
