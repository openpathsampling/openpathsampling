'''
@author David W.H. Swenson
'''
import os, os.path

from xyztranslator import XYZTranslator
from TrajFile import TrajFile, TrajFrame, trajs_equal
from nose.tools import assert_equal
from nose.tools import assert_not_equal
from nose.tools import raises

class testXYZTranslator(object):
    '''Technically, this is more about integration testing than unit testing
    per se, but the point is that these tests should add up to us being able
    to round-trip between the formats.'''

    def setUp(self):
        self.translator=XYZTranslator()
        # fill translator's trajfile object with fake data
        labels = ["A", "B"]
        mass = [1.5, 1.8]

        pos1 = [ [0.5,0.5,0.5], [0.6,0.6,0.6] ]
        vel1 = [ [0.5,0.5,0.5], [0.6,0.6,0.6] ]
        f1 = TrajFrame()
        f1.natoms = 2
        f1.mass = mass
        f1.labels = labels
        f1.pos = pos1
        f1.vel = vel1

        pos2 = [ [0.7,0.7,0.7], [0.8,0.8,0.8] ]
        vel2 = [ [0.7,0.7,0.7], [0.8,0.8,0.8] ]
        f2 = TrajFrame()
        f2.natoms = 2
        f2.mass = mass
        f2.labels = labels
        f2.pos = pos2
        f2.vel = vel2

        self.translator.trajfile = TrajFile()
        self.translator.trajfile.frames = [ f1, f2 ]


    def teardown(self):
        if os.path.isfile("test.nc"): os.remove("test.nc")
        del self.translator

    def test_guess_fname_format(self):
        '''
        Obtain the filename format for XYZ trajfiles from a given filename.
        '''
        test1="mytraj0003.xyz"
        result1="mytraj%04d.xyz"
        test2="mytraj00_0003.xyz"
        result2="mytraj00_%04d.xyz"
        test3="mytraj00_0003acc.xyz"
        result3="mytraj00_%04dacc.xyz"
        fail3="mytraj%02d_0003acc.xyz"
        assert_equal(self.translator.guess_fname_format(test1), result1)
        assert_equal(self.translator.guess_fname_format(test2), result2)
        assert_equal(self.translator.guess_fname_format(test3), result3)
        assert_not_equal(self.translator.guess_fname_format(test3),fail3)

    def test_trajfile_trajectory_trajfile(self):
        '''Integration test: round trip TrajFile->trajectory->TrajFile'''
        self.translator.trajfile2trajectory(self.translator.trajfile)
        # TODO assert that we have the appropriate number of trajectories
        # and the appropriate number of frames
        oldtrajfile = self.translator.trajfile
        self.translator.trajfile = None
        self.translator.trajectory2trajfile(self.translator.trajectory)
        # the two TrajFile objects shouldn't be the same object
        assert_not_equal(self.translator.trajfile, oldtrajfile)
        # but they should be equal on a deep content comparison
        assert_equal(trajs_equal(self.translator.trajfile, oldtrajfile),True)

