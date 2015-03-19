import openpathsampling as paths
from paths.dynamics_engine import DynamicsEngine
from paths.todict import restores_as_full_object

import subprocess

import logging

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
    }


    def __init__(self, options, template):
        pass

    @property
    def current_snapshot(self):
        pass

    @current_snapshot.setter
    def current_snapshot(self, snap):
        pass

    def generate_next_frame(self):
        pass

    def read_frame_from_file(self, frame_num):
        pass

    def start(self, snapshot=None):
        super(ExternalEngine, self).start(snapshot)
        pass

    def stop(self, trajectory):
        super(ExternalEngine, self).stop(trajectory)
        pass
