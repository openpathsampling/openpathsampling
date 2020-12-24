# UUID recognition and encoding #####################################
# Things in here might be modified for performance optimization. In
# particular, it might be worth using a string representation of the UUID
# whenever possible (dicts with string keys have a special fast-path)

import re
import sys
if sys.version_info > (3, ):
    unicode = str
    long = int

def has_uuid(obj):
    return hasattr(obj, '__uuid__')


def get_uuid(obj):
    # TODO: I can come up with a better string encoding than this
    try:
        return str(obj.__uuid__)
    except AttributeError as e:
        if obj is None:
            return obj
        else:
            raise e

def set_uuid(obj, uuid):
    obj.__uuid__ = long(uuid)


def encode_uuid(uuid):
    return "UUID(" + str(uuid) + ")"


def decode_uuid(uuid_str):
    return uuid_str[5:-1]


# use the regular expression when looking through an entire JSON string; use
# the is_uuid_string method for individual objects
encoded_uuid_re = re.compile("UUID\((?P<uuid>[\-]?[0-9]+)\)")


def is_uuid_string(obj):
    return (
        isinstance(obj, (str, unicode))
        and obj[:5] == 'UUID(' and obj[-1] == ')'
    )


