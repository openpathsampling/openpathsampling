from .serialization import *
import pytest

from .serialization_helpers import get_uuid
from .test_utils import (LoadingStorageMock, all_objects, toy_uuid_maker,
                         MockUUIDObject)

class TestGenericLazyLoader(object):
    def setup(self):
        original_and_class = {
            'normal': (all_objects['int'], MockUUIDObject)
            # TODO: add iterable and mappable classes
        }
        self.originals = {k: v[0] for k, v in original_and_class.items()}
        original_class = {k: v[1] for k, v in original_and_class.items()}
        uuid_dict = {get_uuid(obj): obj for obj in self.originals.values()}
        self.storage = LoadingStorageMock(uuid_dict)
        self.proxies = {
            proxy_type: GenericLazyLoader(
                get_uuid(self.originals[proxy_type]),
                original_class[proxy_type],
                self.storage
            )
            for proxy_type in self.originals.keys()
        }

    @pytest.mark.parametrize('proxy_type',
                             ['normal'])#, 'iterable', 'mappable'])
    def test_init(self, proxy_type):
        proxy = self.proxies[proxy_type]
        original = self.originals[proxy_type]
        assert proxy._loaded_object is None
        assert get_uuid(proxy) == get_uuid(original)
        assert proxy.storage == self.storage
        assert isinstance(original, proxy.class_)
        assert proxy._loaded_object is None

    @pytest.mark.parametrize('proxy_type',
                             ['normal'])#, 'iterable', 'mappable'])
    def test_load(self, proxy_type):
        proxy = self.proxies[proxy_type]
        original = self.originals[proxy_type]
        assert proxy._loaded_object is None
        proxy.load()
        assert proxy._loaded_object == original

        # now we check that a second load() doesn't require storage
        proxy.storage = LoadingStorageMock({})
        try:
            proxy.load()
        except RuntimeError:
            raise AssertionError("Proxy required revisiting storage "
                                 "after load")

    def test_load_error(self):
        good_proxy = self.proxies['normal']
        storage = LoadingStorageMock({})
        bad_proxy = GenericLazyLoader(get_uuid(good_proxy),
                                      good_proxy.class_,
                                      storage)
        assert bad_proxy._loaded_object is None
        # TODO: current mock storage raises KeyError instead of returning
        # None; check the behavior of the actual storage -- maybe this isn't
        # possible and the relevant lines need to be removed?
        # with pytest.raises(RuntimeError):
            # bad_proxy.load()
