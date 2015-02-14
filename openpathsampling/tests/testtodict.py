'''
@author: Jan-Hendrik Prinz
'''

import os
import numpy as np

import os
from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import (true_func, data_filename,
                          assert_equal_array_array,
                          assert_not_equal_array_array)


import openpathsampling as paths

from openpathsampling.todict import ObjectJSON, class_list

class testSampleSet(object):
    def setup(self):
        self.json = ObjectJSON()
        self.filename = data_filename("toy_tis.nc")
        self.storage = paths.storage.Storage(filename=self.filename, mode='a')

    def test_todict_full(self):
        # this uses a storage. Get all contained objects, turn them from json
        # to an object and back

        for store_name in self.storage.list_stores():
#            print '%s' % store_name
            store = getattr(self.storage, store_name)
            if hasattr(store.content_class, 'creatable'):
                for obj in store:
                    print 'store %s object %d' % (store_name, obj.idx[self.storage])
                    json1 = obj.json
                    obj_clone = self.json.build(obj.json)
                    json2 = self.json.simplify(obj_clone)
                    assert(json1 == json2)


    def test_todict_stub(self):
        # this uses a storage. Get all contained objects, turn them from json
        # to an object and back

        for store_name in self.storage.list_stores():
#            print '%s' % store_name
            store = getattr(self.storage, store_name)
            if hasattr(store.content_class, 'dictable'):
                for obj in store:
                    print 'store %s object %d' % (store_name, obj.idx[self.storage])
                    json1 = obj.json
                    obj_clone = self.json.build(obj.json)
                    json2 = self.json.simplify(obj_clone)
                    assert(json1 == json2)

