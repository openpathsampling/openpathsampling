from nose.tools import assert_equal, assert_not_equal, assert_is, raises
from nose.plugins.skip import Skip, SkipTest
from test_helpers import CallIdentity, raises_with_message_like

import unittest

import openpathsampling.treelogic as tree
from openpathsampling.treelogic import TupleTree as tup


def setUp():
    global op_id, volA, volB, volC, volD, volA2

class testTupleTree(object):
    def setUp(self):
        # setup a tree like
        # ('a', ('b', ('a', ('d')), ('e', ('d'))))
        self.t = tup(['a', tup(['b', tup(['a', tup(['d'])]), tup(['e', tup(['d'])])])])
        self.t2 = tup(['a', tup(['b', tup(['a', tup(['d'])])])])


    def test_treeprint(self):
        structure = self.t.treeprint().replace('\n', '')
        assert(structure == 'a +- b     +- a     |   +- d     +- e         +- d')