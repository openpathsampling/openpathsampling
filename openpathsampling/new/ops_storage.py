import storage
import sql_backend

from serialization_helpers import to_json_obj as serialize_sim
from serialization_helpers import from_json_obj as deserialize_sim

import openpathsampling as paths

from storage import ClassInfo

ops_schema = {
    'samples': [('trajectory', 'lazy'), ('ensemble', 'uuid'),
                ('replica', 'int'),
                # in my opinion, the next 3 should be removed
                ('parent', 'lazy'), ('bias', 'float'),
                ('mover', 'uuid')],
    'sample_sets': [('samples', 'list_uuid'), ('movepath', 'lazy')],
    'trajectories': [('snapshots', 'list_uuid')],
    'move_changes': [('mover', 'uuid'), ('details', 'lazy'), ('cls', 'str'),
                     ('subchanges', 'list_uuid'), ('samples', 'list_uuid'),
                     ('input_samples', 'list_uuid')],
    'steps': [('change', 'uuid'), ('active', 'uuid'), ('previous', 'lazy'),
              ('simulation', 'uuid'), ('mccycle', 'int')],
    'details': [('json', 'json')],
    'simulation_objects': [('json', 'json_obj'), ('class_idx', 'int')]
}

ops_schema_sql_metadata = {}

class OPSClassInfoContact(storage.ClassInfoContainer):
    def is_special(self, item):
        return isinstance(item, paths.BaseSnapshot)

    def special_lookup(self, item):
        if isinstance(item, paths.BaseSnapshot):
            return (item.engine, item.__class__)

ops_class_info = OPSClassInfoContainer(
    default_info=ClassInfo('simulation_objects', cls=StorableObject,
                           serializer=serialize_sim,
                           deserializer=deserialize_sim),
    class_info_list=[
        ClassInfo(table='samples', cls=paths.Sample)
        ClassInfo(table='sample_sets', cls=paths.SampleSet)
        ClassInfo(table='trajectories', cls=paths.Trajectory)
        ClassInfo(table='move_changes', cls=paths.MoveChange)
        ClassInfo(table='steps', cls=paths.MCStep)
        ClassInfo(table='details', cls=paths.Details)
    ]
)


class OPSStorage(storage.GeneralStorage):
    def __init__(self, backend, schema, class_info, fallbacks=None):
        super(OPSStorage, self).__init__(backend, schema, class_info,
                                         fallbacks)
        self.n_snapshot_types = 0

