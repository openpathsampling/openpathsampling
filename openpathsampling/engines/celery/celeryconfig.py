# Register your new serializer methods into kombu
from kombu.serialization import register

import openpathsampling.netcdfplus.dictify as dfy

_ops_to = dfy.CachedUUIDObjectJSON()
# _ops_from = dfy.CachedUUIDObjectJSON()

register('opsjson', _ops_to.to_json, _ops_to.from_json,
    content_type='application/x-opsjson',
    content_encoding='utf-8')

# Tell celery to use your new serializer:
CELERY_ACCEPT_CONTENT = ['opsjson']
CELERY_TASK_SERIALIZER = 'opsjson'
CELERY_RESULT_SERIALIZER = 'opsjson'
