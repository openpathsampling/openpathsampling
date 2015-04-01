import openpathsampling as paths
from openpathsampling import DynamicsEngine
from openpathsampling import restores_as_full_object
import os

import logging

import psutil
import signal
import shlex
import time

@restores_as_full_object
class ExternalEngine(DynamicsEngine):
    """
    Generic object to handle arbitrary external engines. 

    Typically, this will be subclassed for any given engine. As written, it
    will work with the trivial `engine.c` developed for testing purposes.
    """

    default_options = {
        'n_frames_max' : 10000,
        'name_prefix' : "test",
        'default_sleep_ms' : 10,
        'engine_sleep' : 100
    }

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
        self.killsig = signal.SIGKILL
        self.traj_num = -1

    @property
    def current_snapshot(self):
        return self._current_snapshot

    @current_snapshot.setter
    def current_snapshot(self, snap):
        self._current_snapshot = snap
        # for a real engine, this should also write frame to whatever
        # file is input for the trajectory

    def generate_next_frame(self):
        # should be completely general
        next_frame_found = False
        while not next_frame_found:
            try:
                next_frame = self.read_frame_from_file(self.frame_num)
                next_frame_found = True
            except: #TODO what error to throw here?
                time.sleep(self.sleep_ms/1000.0)
                # TODO: optimize sleep time
        self.current_snapshot = next_frame
        return self.current_snapshot

    def start(self, snapshot=None):
        super(ExternalEngine, self).start(snapshot)
        self.traj_num += 1
        self.set_input_filename(self.traj_num)
        self.set_output_filename(self.traj_num)
        if snapshot is not None:
            self.current_snapshot = snapshot
        self.write_frame_to_file(self.input_file, self.current_snapshot)

        cmd = shlex.split(self.engine_command())
        try:
            self.proc = psutil.Popen(shlex.split(self.engine_command()),
                                         preexec_fn=os.setsid
                                        )
        except OSError:
            pass #TODO: need to handle this, but do what?

    def stop(self, trajectory):
        super(ExternalEngine, self).stop(trajectory)
        self.proc.send_signal(self.killsig)
        self.proc.wait() # wait for the zombie to die

    # FROM HERE ARE THE FUNCTIONS TO OVERRIDE IN SUBCLASSES:
    def read_frame_from_file(self, filename, frame_num):
        """Reads given frame number from file, and returns snapshot
        """
        # main part that needs to be rewritten for each engine
        pass

    def write_frame_to_file(self, filename, snapshot):
        """Writes given snapshot to appropriate input file.
        """
        # also needs to be rewritten for each engine
        pass

    def set_input_filename(self, number):
        """Sets the filename for the input to trajectory number `number`
        """
        self.input_file = self.name_prefix + str(number) + ".inp"

    def set_output_filename(self, number):
        """Sets the filename for the output from trajectory number `number`
        """
        self.output_file = self.name_prefix + str(number) + ".out"

    def engine_command(self):
        """Generates a string for the command to run the engine."""
        return "engine " + str(self.engine_sleep) + " " + str(self.output_file)


