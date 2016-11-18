import json
from datetime import datetime
from time import mktime

import openpathsampling.netcdfplus as npl

class MyEncoder(json.JSONEncoder):
    def __init__(self, *args, **kwargs):
        super(MyEncoder, self).__init__(*args, **kwargs)
        self._ops_dictifier = npl.ObjectJSON()

    def default(self, obj):
        try:
            return self._ops_dictifier.simplify(obj)
        except:
        # if isinstance(obj, npl.StorableObject):
        #     return {
        #         '__type__': '__datetime__',
        #         'epoch': int(mktime(obj.timetuple()))
        #     }
        # else:
            return json.JSONEncoder.default(self, obj)

def my_decoder(obj):
    # if '__type__' in obj:
    #     if obj['__type__'] == '__datetime__':
    #         return datetime.fromtimestamp(obj['epoch'])

    obj = self._ops_dictifier.build(obj)

    return obj

# Encoder function
def my_dumps(obj):
    return json.dumps(obj, cls=MyEncoder)

# Decoder function
def my_loads(obj):
    return json.loads(obj, object_hook=my_decoder)