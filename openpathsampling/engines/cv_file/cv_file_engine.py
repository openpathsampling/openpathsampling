import math
import pandas as pd
import numpy as np
from openpathsampling.engines.external_engine import ExternalEngine
from openpathsampling.engines.toy import ToySnapshot
import openpathsampling as paths

# inspired by a question Gerhard Hummer asked at a Lorentz Center / ECAM
# state of the art meeting.

class CVFileEngine(ExternalEngine):

    _default_options = {
        'path_to_engine': "a.out",
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

        self.velocities = np.zeros((n_atoms, n_spatial))  # HACK

    def trajectory_from_file(self, filename):
        df = pd.read_table(filename, delim_whitespace=True)
        snaplist = [ToySnapshot(coordinates=df.iloc[i].values,
                                velocities=self.velocities)
                    for i in df.index]
        return paths.Trajectory(snaplist)

    def read_frame_from_file(self, filename, frame_num):
        # not the most efficient, but a quick hack
        df = pd.read_table(filename, delim_whitespace=True)
        return ToySnapshot(coordinates=df.iloc[df.index[frame_num]],
                           velocities=self.velocities)

    def write_frame_to_file(self, filename, snapshot, mode="a"):
        flattened = sum(snapshot.coordinates.to_list(), [])
        df = pd.DataFrame(data=flattened)
        df.to_csv(filename, sep=" ", header=False, index=False, mode=mode)

    def engine_command(self):
        return (self.path_to_engine + " " + str(self.input_file) + " > " +
                self.output_file)
