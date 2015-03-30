import openpathsampling as paths
from openpathsampling import DynamicsEngine
from openpathsampling import restores_as_full_object

import logging

import subprocess
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
        'engine_command' : "engine 100 engine.out",
        'trajectory_file' : "engine.out"
    }

    def __init__(self, options, template):
        # needs to be overridden for each engine
        options = {
            'n_spatial' : 1,
            'n_atoms' : 1
        }
        super(ExternalEngine, self).__init__(options=options)
        self.template = template

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
        pass

    def read_frame_from_file(self, frame_num):
        # main part that needs to be rewritten for each engine
        pass

    def set_input_file(self, fname):
        self.input_file = fname

    def set_output_file(self, fname):
        self.output_file = fname

    def engine_command(self):
        return "engine " + str(delay_time) + " " + str(self.output_file)

    def start(self, snapshot=None):
        super(ExternalEngine, self).start(snapshot)
        try:
            self.proc = subprocess.Popen(shlex.split(self.engine_command()))
        except OSError:
            pass

    def stop(self, trajectory):
        super(ExternalEngine, self).stop(trajectory)
        self.proc.kill()
