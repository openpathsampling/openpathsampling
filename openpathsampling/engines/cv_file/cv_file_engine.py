import math
import pandas as pd
from openpathsampling.engines.external_engine import ExternalEngine
import openpathsampling as paths

# inspired by a question Gerhard Hummer asked at a Lorentz Center / ECAM
# state of the art meeting.

class CVFileEngine(ExternalEngine):

    _default_options = {
        'engine_command': "a.out",
        'cv_names': []
    }
    def __init__(self, options=None):
        # FUTURE: add `format` (csv, ssv, etc)
        try:
            n_atoms = options['n_atoms']
        except KeyError:
            n_atoms = 1
        try:
            cv_names = options['cv_names']
        except KeyError:
            cv_names = self._default_options['cv_names']
        try:
            n_spatial = options['n_spatial']
        except KeyError:
            if len(cv_names) == 0:
                n_spatial = 1
            else:
                n_spatial = math.ceil(float(len(cv_names)) / n_atoms)

        template = np.array([[0.0]*len(cv_name_list)])
        if options is None:
            options = {}
        super(CVFileEngine, self).__init__(options=options, template)
        self.cv = {}
        for i in range(len(cv_names)):
            name = cv_names[i]
            atom_i = i / n_atoms
            spatial_i = i % n_spatial
            self.cv[name] = paths.FunctionCV(name, lambda s: s.xyz[i])

    def trajectory_from_file(self, filename):
        pass

    def read_frame_from_file(self, filename, frame_num):
        pass

    def write_frame_to_file(self, filename, snapshot, mode="a"):
        pass

    def engine_command(self):
        pass
