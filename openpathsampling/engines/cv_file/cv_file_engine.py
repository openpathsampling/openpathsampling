import math
import pandas as pd
import numpy as np
from openpathsampling.engines.external_engine import ExternalEngine
from openpathsampling.engines.toy import ToySnapshot
import openpathsampling as paths

# inspired by a question Gerhard Hummer asked at a Lorentz Center / E-CAM
# state of the art meeting.

class CVFileEngine(ExternalEngine):

    _default_options = {
        'path_to_engine': "a.out",
        'cv_names': ['x'],
        'n_frames_max' : 10000,
        'name_prefix' : "test",
        'default_sleep_ms' : 100,
        'auto_optimize_sleep' : True,
        'engine_sleep' : 100,
        'engine_directory' : "",
        'n_spatial' : 1,
        'n_atoms' : 1
    }
    def __init__(self, options=None):

        # FUTURE: add `format` (csv, ssv, etc)
        if options is None:
            options = {}
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
                n_spatial = int(math.ceil(float(len(cv_names)) / n_atoms))

        options['n_atoms'] = n_atoms
        options['n_spatial'] = n_spatial

        template = np.array([[0.0]*len(cv_names)])
        super(CVFileEngine, self).__init__(options=options,
                                           template=template)

        # explanation for the strange i=i is because lambdas are annoying
        # when defined in loops: http://stackoverflow.com/questions/938429
        cvs = [paths.FunctionCV(self.cv_names[i],
                                lambda s, i=i: s.xyz.flatten()[i])
               for i in range(len(self.cv_names))]
        self.cv = {cv.name: cv for cv in cvs}

        self.shape = (n_atoms, n_spatial)
        self.velocities = np.zeros(self.shape)  # HACK

    def trajectory_from_file(self, filename):
        df = pd.read_table(filename, delim_whitespace=True, header=None)
        snaplist = [
            ToySnapshot(
                coordinates=df.iloc[i].values.reshape(self.shape),
                velocities=self.velocities
            )
            for i in df.index
        ]
        return paths.Trajectory(snaplist)

    def read_frame_from_file(self, filename, frame_num):
        # not the most efficient, but a quick hack
        df = pd.read_table(filename, delim_whitespace=True, header=None)
        coords = df.iloc[df.index[frame_num]].values.reshape(self.shape)
        return ToySnapshot(coordinates=coords,
                           velocities=self.velocities)

    def write_frame_to_file(self, filename, snapshot, mode="a"):
        flattened = sum(snapshot.coordinates.to_list(), [])
        df = pd.DataFrame(data=flattened)
        df.to_csv(filename, sep=" ", header=False, index=False, mode=mode)

    def engine_command(self):
        # TODO: self.input_file doesn't exist
        return (self.path_to_engine + " " + str(self.input_file) + " > " +
                self.output_file)
