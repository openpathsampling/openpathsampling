import re
import numpy as np
ndarray_re = \
        re.compile("ndarray\.(?P<dtype>[a-z0-9]+)(?P<shape>\([0-9\,\ ]+\))")

def parse_ndarray_type(type_name):
    m_ndarray = ndarray_re.match(type_name)
    if m_ndarray:
        dtype = getattr(np, m_ndarray.group('dtype'))
        shape = tuple(map(int, m_ndarray.group('shape')[1:-1].split(',')))
        return dtype, shape
    return None


uuid_types = ['uuid', 'lazy']
list_uuid = ['list_uuid']
builtin_types = ['str', 'int', 'float']
ndarray_types = ['ndarray']

all_types = uuid_types + list_uuid + builtin_types + ndarray_types


