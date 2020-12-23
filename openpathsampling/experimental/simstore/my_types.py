import re
import numpy as np
ndarray_re = re.compile(
    "ndarray\.(?P<dtype>[a-z0-9]+)(?P<shape>\([0-9\,\ ]+\))"
)

def parse_ndarray_type(type_name):
    m_ndarray = ndarray_re.match(type_name)
    if m_ndarray:
        dtype = getattr(np, m_ndarray.group('dtype'))
        shape = tuple(map(int, m_ndarray.group('shape')[1:-1].split(',')))
        return dtype, shape
    return None

# TODO: this needs to be set up in a way to make it extensible (without
# editing core code)
def backend_registration_type(type_name):
    backend_type = type_name
    ndarray_info = parse_ndarray_type(type_name)
    if parse_ndarray_type(type_name):
        backend_type = 'ndarray'
    return backend_type



uuid_types = ['uuid', 'lazy']
uuid_list_types = ['list_uuid']
builtin_types = ['str', 'int', 'float']
ndarray_types = ['ndarray']

all_types = uuid_types + uuid_list_types + builtin_types + ndarray_types
all_uuid_types = ['uuid', 'lazy', 'list_uuid']
json_obj_types = ['json', 'json_obj']


