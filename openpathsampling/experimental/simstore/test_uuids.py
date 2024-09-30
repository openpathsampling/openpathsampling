import pytest

from .uuids import *
from .test_utils import create_test_objects, MockUUIDObject

all_objects = create_test_objects()

@pytest.mark.parametrize('obj', list(all_objects.values()))
def test_has_uuid(obj):
    assert has_uuid(obj)

def test_has_uuid_no_uuid():
    assert not has_uuid(10)
    assert not has_uuid('foo')

@pytest.mark.parametrize('name,obj', list(all_objects.items()))
def test_get_uuid(name, obj):
    assert get_uuid(obj) == str(int(hash(name)))

def test_get_uuid_none():
    assert get_uuid(None) is None

def test_get_uuid_error():
    with pytest.raises(AttributeError):
        get_uuid(10)

def test_set_uuid():
    obj = MockUUIDObject(name="test", normal_attr=10)
    fake_uuid = 100
    assert get_uuid(obj) != str(fake_uuid)
    set_uuid(obj, fake_uuid)
    assert get_uuid(obj) == str(fake_uuid)

def test_encode_uuid():
    obj = MockUUIDObject(name="test", normal_attr=10)
    set_uuid(obj, 100)
    assert encode_uuid("UUID(100)")


