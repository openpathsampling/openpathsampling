from .mdtraj_json import *
import pytest

import numpy as np
import numpy.testing as npt

from ..simstore.custom_json import bytes_codec, numpy_codec, custom_json_factory
from ..simstore.test_custom_json import CustomJSONCodingTest

from openpathsampling.tests.test_helpers import data_filename

class MDTrajCodingTest(CustomJSONCodingTest):
    def setup(self):
        if not HAS_MDTRAJ:
            pytest.skip()

        self.filename = data_filename('ala_small_traj.pdb')

    def test_default(self):
        # custom for handling numpy
        for (obj, dct) in zip(self.objs, self.dcts):
            default = self.codec.default(obj)
            numpy_attrs = [attr for attr, val in dct.items()
                           if isinstance(val, np.ndarray)]
            other_attrs = [attr for attr, val in dct.items()
                           if not isinstance(val, np.ndarray)]
            for attr in numpy_attrs:
                npt.assert_array_equal(default[attr], dct[attr])
            for attr in other_attrs:
                assert default[attr] == dct[attr]

    def test_round_trip(self):
        codecs = [numpy_codec, bytes_codec] + mdtraj_codecs
        encoder, decoder = custom_json_factory(codecs)
        self._test_round_trip(encoder, decoder)


class TestTopologyCoding(MDTrajCodingTest):
    def setup(self):
        super(TestTopologyCoding, self).setup()
        self.codec = top_codec
        top = md.load(self.filename).topology
        dataframe, bonds = top.to_dataframe()
        self.objs = [top]
        self.dcts = [{
            '__class__': 'Topology',
            '__module__': 'mdtraj.core.topology',
            'atoms': dataframe.to_json(),
            'bonds': bonds
        }]


class TestTrajectoryCoding(MDTrajCodingTest):
    def setup(self):
        super(TestTrajectoryCoding, self).setup()
        self.codec = traj_codec
        traj = md.load(self.filename)
        self.objs = [traj]
        self.dcts = [{
            '__class__': 'Trajectory',
            '__module__': 'mdtraj.core.trajectory',
            'xyz': traj.xyz,
            'topology': traj.topology,
            'time': traj.time,
            'unitcell_lengths': traj.unitcell_lengths,
            'unitcell_angles': traj.unitcell_angles
        }]



