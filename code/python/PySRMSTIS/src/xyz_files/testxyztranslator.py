'''
@author David W.H. Swenson
'''
import os
import shutil

from xyztranslator import XYZTranslator
from TrajFile import TrajFile, TrajFrame, trajs_equal
from nose.tools import assert_equal
from nose.tools import assert_not_equal
from nose.tools import raises

class testXYZTranslator(object):
    '''Technically, this is more about integration testing than unit testing
    per se, but the point is that these tests should add up to us being able
    to round-trip between the formats.
    '''

    def setUp(self):
        self.translator=XYZTranslator()
        # fill translator's trajfile object with fake data
        labels = ["A", "B", "A"]
        mass = [1.5, 1.8, 1.5]

        pos1 = [ [0.5,0.5,0.5], [0.6,0.6,0.6], [0.4,0.4,0.4] ]
        vel1 = [ [1.5,1.5,1.5], [1.6,1.6,1.6], [1.4,1.4,1.4] ]
        f1 = TrajFrame()
        f1.natoms = 3
        f1.mass = mass
        f1.labels = labels
        f1.pos = pos1
        f1.vel = vel1

        pos2 = [ [0.7,0.7,0.7], [0.8,0.8,0.8], [0.9,0.9,0.9] ]
        vel2 = [ [1.7,1.7,1.7], [1.8,1.8,1.8], [1.9,1.9,1.9] ]
        f2 = TrajFrame()
        f2.natoms = 3
        f2.mass = mass
        f2.labels = labels
        f2.pos = pos2
        f2.vel = vel2

        self.translator.trajfile = TrajFile()
        self.translator.trajfile.frames = [ f1, f2 ]


    def teardown(self):
        if os.path.isfile("test.nc"): os.remove("test.nc")
        if os.path.isdir("mydir00"): shutil.rmtree("mydir00")
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
        self.translator.outfile = ("test.nc", "nc")
        self.translator.trajectory = self.translator.trajfile2trajectory(self.translator.trajfile)
        # TODO assert that we have the appropriate number of trajectories
        # and the appropriate number of frames
        oldtrajfile = self.translator.trajfile
        self.translator.trajfile = None
        self.translator.trajfile = self.translator.trajectory2trajfile(self.translator.trajectory)
        # the two TrajFile objects shouldn't be the same object
        assert_not_equal(self.translator.trajfile, oldtrajfile)
        # but they should be equal on a deep content comparison
        assert_equal(trajs_equal(self.translator.trajfile, oldtrajfile),True)

    def test_set_infile_outfile(self):
        # get filenames for xyz->nc
        testin=["sometest_00_000003.xyz"]
        testout="sometest_00_000003.nc"
        self.translator.set_infile_outfile(testin, testout)
        assert_equal(self.translator.infiles, testin)
        assert_equal(self.translator.outfile, testout)
        assert_equal(self.translator.intype, "xyz")
        assert_equal(self.translator.outtype, "nc")
        # get filenames for nc->xyz
        testin=["sometest_00_000003.nc"]
        testout="sometest_00_000003.xyz"
        self.translator.set_infile_outfile(testin, testout)
        assert_equal(self.translator.infiles, testin)
        assert_equal(self.translator.outfile, testout)
        assert_equal(self.translator.intype, "nc")
        assert_equal(self.translator.outtype, "xyz")
        # split off directory for nc->xyz dir
        testin=["sometest.nc"]
        testout="mydir00/myfile_00_%06d.xyz"
        self.translator.set_infile_outfile(testin, testout)
        assert_equal(self.translator.infiles, testin)
        assert_equal(self.translator.outfile, testout)
        assert_equal(self.translator.intype, "nc")
        assert_equal(self.translator.outtype, "xyz")
        assert_equal(os.path.isdir("mydir00"), True)
        # get multiple files for xyz dir->nc
        testin=["mydir00/myfile_00_000000.xyz", "mydir00/myfile_00_000001.xyz"]
        testout="sometest.nc"
