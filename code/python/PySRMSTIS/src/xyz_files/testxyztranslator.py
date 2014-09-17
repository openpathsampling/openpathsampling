'''
@author David W.H. Swenson
'''

from xyztranslator import XYZTranslator
from TrajFile import trajs_equal
from nose.tools import assert_equal
from nose.tools import assert_not_equal
from nose.tools import raises

class testXYZTranslator(object):
    '''Technically, this is more about integration testing than unit testing
    per se, but the point is that these tests should add up to us being able
    to round-trip between the formats.'''

    @classmethod
    def setup_class(myclass):
        pass

    @classmethod
    def teardown_class(myclass):
        pass

    def setUp(self):
        # TODO: this is where you add fake data
        self.translator=XYZTranslator()
        pass

    def teardown(self):
        pass

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
        assert_equal(self.translator.guess_fname_format(test1), result1)
        assert_equal(self.translator.guess_fname_format(test2), result2)
        assert_equal(self.translator.guess_fname_format(test3), result3)
        assert_not_equal(self.translator.guess_fname_format(test3, 
                        "mytraj%02d_0003acc.xyz"))

    def test_trajfile_storage_trajfile(self):
        '''Integration test: '''
        self.translator.trajfile2trajectory(self.translator.trajfile)
        oldtrajfile = self.translator.trajfile
        self.translator.trajfile = None
        self.translator.trajector2trajfile(self.translator.trajectory)
        # the two TrajFile objects shouldn't be the same object
        assert_not_equal(self.translator.trajfile, oldtrajfile)
        # but they should be equal on a deep content comparison
        assert_equal(trajs_equal(self.translator.trajfile, oldtrajfile))

