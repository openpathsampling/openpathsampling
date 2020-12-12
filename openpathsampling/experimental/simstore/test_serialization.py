from .serialization import *
from .proxy import *  # TODO: move this elsewhere
import pytest

from .serialization_helpers import get_uuid, default_find_uuids
from .test_utils import (LoadingStorageMock, all_objects, toy_uuid_maker,
                         MockUUIDObject, MockSimulationObject, MockBackend)

from . import class_info

from .attribute_handlers import DEFAULT_HANDLERS


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
        pytest.skip()
        # TODO: current mock storage raises KeyError instead of returning
        # None; check the behavior of the actual storage -- maybe this isn't
        # possible and the relevant lines need to be removed?
        # with pytest.raises(RuntimeError):
            # bad_proxy.load()

    def test_getattr(self):
        proxy = self.proxies['normal']
        original = self.originals['normal']
        assert proxy._loaded_object is None
        assert proxy.normal_attr == original.normal_attr
        assert proxy._loaded_object == original

    def test_serialize_proxy(self):
        proxy = self.proxies['normal']
        original = self.originals['normal']
        dct = proxy.to_dict()
        obj = proxy.class_.from_dict(dct)
        assert obj.__uuid__ == proxy.__uuid__ == original.__uuid__
        assert type(obj) == type(original)
        assert obj.normal_attr == proxy.normal_attr == original.normal_attr

    def test_save_proxy(self):
        # TODO: storing a lazy proxy to storage should actually store it as
        # if it the original object (i.e., store to DB correctly ) -- note
        # that only schema-based storage can invoke the lazy proxies, so we
        # this shouldn't risk saving a simulation object
        pytest.skip()


class TestProxyObjectFactory(object):
    def setup(self):
        self.storage = LoadingStorageMock({get_uuid(obj): obj
                                           for obj in all_objects.values()})

    def test_make_lazy(self):
        factory = ProxyObjectFactory(self.storage, None)
        assert len(factory.lazy_classes) == 0
        original = all_objects['int']
        uuid = get_uuid(original)
        proxy = factory.make_lazy(MockUUIDObject, uuid)
        assert len(factory.lazy_classes) == 1
        assert list(factory.lazy_classes.keys()) == [MockUUIDObject]
        assert isinstance(proxy, MockUUIDObject)
        assert proxy._loaded_object is None
        assert proxy.normal_attr == original.normal_attr

    def test_make_all_lazies(self):
        backend = MockBackend()
        obj = all_objects['obj']
        uuid = get_uuid(obj)
        lazy_rows = backend.load_uuids_table([uuid])
        lazies = {'mock': lazy_rows}

        # this includes integration with serialization schema
        schema = {'mock': [('obj_attr', 'lazy')]}
        sim_info = class_info.ClassInfo('simulation_objects',
                                        cls=MockSimulationObject,
                                        find_uuids=default_find_uuids)
        mock_info = class_info.ClassInfo(table='mock', cls=MockUUIDObject)
        mock_info.set_defaults(schema, DEFAULT_HANDLERS)
        serialization_schema = class_info.SerializationSchema(
            default_info=sim_info,
            schema=schema,
            class_info_list=[mock_info]
        )
        factory = ProxyObjectFactory(self.storage, serialization_schema)
        lazy_objs = factory.make_all_lazies(lazies)
        assert len(lazy_objs) == 1
        proxy = lazy_objs[uuid]
        assert isinstance(proxy, MockUUIDObject)
        assert proxy._loaded_object is None
        assert proxy.obj_attr.normal_attr == 5
        assert proxy._loaded_object is not None
