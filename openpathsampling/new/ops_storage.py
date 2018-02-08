import storage
import sql_backend

from serialization_helpers import to_json_obj as serialize_sim
from serialization_helpers import from_json_obj as deserialize_sim

import openpathsampling as paths
from openpathsampling.netcdfplus import StorableObject

from storage import ClassInfo

import snapshots

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

class OPSClassInfoContainer(storage.ClassInfoContainer):
    def is_special(self, item):
        return isinstance(item, paths.BaseSnapshot)

    def special_lookup_key(self, item):
        if isinstance(item, paths.BaseSnapshot):
            return (item.engine, item.__class__)

ops_class_info = OPSClassInfoContainer(
    default_info=ClassInfo('simulation_objects', cls=StorableObject,
                           serializer=serialize_sim,
                           deserializer=deserialize_sim),
    class_info_list=[
        ClassInfo(table='samples', cls=paths.Sample),
        ClassInfo(table='sample_sets', cls=paths.SampleSet),
        ClassInfo(table='trajectories', cls=paths.Trajectory),
        ClassInfo(table='move_changes', cls=paths.MoveChange),
        ClassInfo(table='steps', cls=paths.MCStep),
        ClassInfo(table='details', cls=paths.Details),
    ]
)


class OPSStorage(storage.GeneralStorage):
    def __init__(self, backend, schema, class_info, fallbacks=None):
        # TODO: this will change to match the current notation
        super(OPSStorage, self).__init__(backend, schema, class_info,
                                         fallbacks)
        self.n_snapshot_types = 0

    @classmethod
    def from_backend(cls, backend, schema=None, class_info=None,
                     fallbacks=None):
        obj = cls.__new__(cls)
        if schema is None:
            scheme = ops_schema
        if class_info is None:
            class_info = ops_class_info
        super(OPSStorage, obj).__init__(backend=backend, schema=schema,
                                        class_info=class_info,
                                        fallbacks=fallbacks)
        obj.n_snapshot_types = 0
        return obj

    def register_from_instance(self, lookup, obj):
        if isinstance(obj, paths.BaseSnapshot):
            schema, class_info_list = snapshots.snapshot_registration_info(
                obj, self.n_snapshot_types
            )
            schema = snapshots.replace_schema_dimensions(
                schema, obj.engine.descriptor
            )
            self.register_schema(schema, class_info_list)
            self.n_snapshot_types += 1

