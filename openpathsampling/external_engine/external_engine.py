import openpathsampling as paths
from openpathsampling import DynamicsEngine
from openpathsampling import restores_as_full_object
import os

import logging

import psutil
import signal
import shlex
import time

import linecache

@restores_as_full_object
class ExternalEngine(DynamicsEngine):
    """
    Generic object to handle arbitrary external engines. 

    Typically, this will be subclassed for any given engine. As written, it
    will work with the trivial `engine.c` developed for testing purposes.
    """

    # TODO: include clever adaptive waiting scheme

    default_options = {
        'n_frames_max' : 10000,
        'name_prefix' : "test",
        'default_sleep_ms' : 10,
        'engine_sleep' : 100
    }

    killsig = signal.SIGTERM

    def __init__(self, options, template):
        # needs to be overridden for each engine
        options = {
            'n_spatial' : 1,
            'n_atoms' : 1
        }
        super(ExternalEngine, self).__init__(options=options)
        self.template = template
        self.sleep_ms = self.default_sleep_ms
        # per engine, you can override this to terminate with the signal of
        # your choice
        self._traj_num = -1

    @property
    def current_snapshot(self):
        return self._current_snapshot

    @current_snapshot.setter
    def current_snapshot(self, snap):
        self._current_snapshot = snap

    def generate_next_frame(self):
        # should be completely general
        next_frame_found = False
        while not next_frame_found:
            next_frame = self.read_frame_from_file(self.frame_num)
            if next_frame == "partial":
                pass # rerun immediately
            elif next_frame is None:
                # TODO: optimize sleep time
                time.sleep(self.sleep_ms/1000.0)
            elif isinstance(next_frame, paths.Snapshot): # success
                self.current_snapshot = next_frame
                next_frame_found = True
                self.frame_num += 1
            else:
                raise RuntimeError("Strange return value from read_next_frame_from_file")
        return self.current_snapshot

    def start(self, snapshot=None):
        super(ExternalEngine, self).start(snapshot)
        self._traj_num += 1
        self.set_filenames(self._traj_num)
        self.write_frame_to_file(self.input_file, self.current_snapshot)

        cmd = shlex.split(self.engine_command())
        try:
            # TODO: add the ability to have handlers for stdin and stdout
            self.proc = psutil.Popen(shlex.split(self.engine_command()),
                                     preexec_fn=os.setsid
                                    )
        except OSError:
            pass #TODO: need to handle this, but do what? Probably reraise

    def stop(self, trajectory):
        super(ExternalEngine, self).stop(trajectory)
        proc = self.who_to_kill()
        proc.send_signal(self.killsig)
        proc.wait() # wait for the zombie to die
        self.cleanup()

    # FROM HERE ARE THE FUNCTIONS TO OVERRIDE IN SUBCLASSES:
    def read_frame_from_file(self, filename, frame_num):
        """Reads given frame number from file, and returns snapshot.
        
        If no frame is available, returns None. If the frame appears to be
        partially written, returns string "partial".
        """
        # in a more complicated case, you'll need to read in all lines from
        # the first associated with frame_num until either (a) you hit EOF
        # or (b) you have a complete frame; then you need to return
        # appropriately 
        first_line = frame_num
        line = linecache.getline(filename, first_line)
        if line is '':
            snap = None
        else:
            try:
                floatval = float(line)
                snap = paths.Snapshot(coordinates=[[floatval]],
                                      velocities=[[1.0]]
                                     )
            except ValueError:
                snap = "partial"
        return snap

    def write_frame_to_file(self, filename, snapshot, mode="a"):
        """Writes given snapshot to file."""
        pass

    def who_to_kill(self):
        """Returns psutil.Process object to send kill signal to.
        
        Might override to send kill signal to a process other than the one
        directly spawned above (e.g., when launching parallel runs)
        """
        # this should only be called if you're about to kill the process; if
        # the process doesn't exist, you shouldn't be killing anything and
        # it will raise an error
        return self.proc

    def cleanup(self):
        """Any cleanup actions to do after the subprocess dies."""
        pass

    def set_filenames(self, number):
        """Sets names for files associated with trajectory `number`"""
        self.input_file = self.name_prefix + str(number) + ".inp" # not used
        self.output_file = self.name_prefix + str(number) + ".out"

    def engine_command(self):
        """Generates a string for the command to run the engine."""
        return "engine " + str(self.engine_sleep) + " " + str(self.output_file)


