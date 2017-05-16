import logging

import mdtraj as md
from openpathsampling.engines import ExternalEngine

class Gromacs5Engine(ExternalEngine):
    _default_options = {
    }

    def read_frame_from_file(self, filename, frame_num):
        pass

    def write_frame_to_file(self, filename, snapshot, mode='a'):
        pass

    def set_filenames(self, number):
        pass

    def engine_command(self):
        pass
    pass
