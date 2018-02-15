import storage
import sql_backend

from serialization_helpers import to_json_obj as serialize_sim
from serialization_helpers import from_json_obj as deserialize_sim
from serialization_helpers import import_class
from serialization_helpers import get_uuid, set_uuid

import openpathsampling as paths
from openpathsampling.netcdfplus import StorableObject

import tools

from serialization import (
    ToDictSerializer, DefaultSerializer, DefaultDeserializer,
    SimulationObjectSerializer
)

from storage import ClassInfo

import snapshots
import logging
logger = logging.getLogger(__name__)

ops_schema = {
    'samples': [('trajectory', 'lazy'), ('ensemble', 'uuid'),
                ('replica', 'int')],
                # in my opinion, the next 3 should be removed
                # ('parent', 'lazy'), ('bias', 'float'),
                # ('mover', 'lazy')],
    # # movepath no longer exists in sample sets?
    'sample_sets': [('samples', 'list_uuid')], #, ('movepath', 'lazy')],
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

class MoveChangeDeserializer(DefaultDeserializer):
    # in general, I think it would be better to reorg MoveChange to only be
    # one class, but this is aimed at fixing problems with reloading
    # MoveChange objects
    def __init__(self, schema, table):
        super(MoveChangeDeserializer, self).__init__(
            schema=schema,
            table=table,
            cls=None
        )

    def __call__(self, uuid, table_dct, cache_list):
        class_name = table_dct.pop('cls')
        cls = import_class('openpathsampling.movechange', class_name)
        dct = self.make_dct(table_dct, cache_list)
        # from here based on JHP's MoveChangeStore._load
        obj = cls.__new__(cls)
        paths.MoveChange.__init__(obj, mover=dct['mover'])
        obj.samples = dct['samples']
        obj.details = dct['details']
        obj.subchanges = dct['subchanges']
        obj.input_samples = dct['input_samples']
        set_uuid(obj, uuid)
        return obj


class OPSClassInfoContainer(storage.ClassInfoContainer):
    special_superclasses = (paths.BaseSnapshot, paths.MoveChange,
                            paths.Details)
    def is_special(self, item):
        return isinstance(item, self.special_superclasses)

    def special_lookup_key(self, item):
        if isinstance(item, paths.BaseSnapshot):
            return (get_uuid(item.engine), item.__class__)
        elif isinstance(item, paths.MoveChange):
            return paths.MoveChange
        elif isinstance(item, paths.Details):
            # TODO: this should be removed, since all Details classes are
            # equivalent -- unfortunately, JHP's LoaderProxy breaks in that
            # case
            return paths.Details

ops_class_info = OPSClassInfoContainer(
    default_info=ClassInfo('simulation_objects', cls=StorableObject,
                           serializer=SimulationObjectSerializer(),
                           deserializer=deserialize_sim),
    class_info_list=[
        ClassInfo(table='samples', cls=paths.Sample),
        ClassInfo(table='sample_sets', cls=paths.SampleSet),
        ClassInfo(table='trajectories', cls=paths.Trajectory),
        ClassInfo(table='move_changes', cls=paths.MoveChange,
                  deserializer=MoveChangeDeserializer(
                      schema=ops_schema,
                      table='move_changes'
                  )),
        ClassInfo(table='steps', cls=paths.MCStep),
        ClassInfo(table='details', cls=paths.Details,
                  serializer=SimulationObjectSerializer(),
                  deserializer=deserialize_sim),
    ]
)

for info in ops_class_info.class_info_list:
    info.set_defaults(ops_schema)


# TODO: add more to these
ops_simulation_classes = {
    'volumes': paths.Volume,
    'ensembles': paths.Ensemble,
    'pathsimulators': paths.PathSimulator,
    'pathmovers': paths.PathMover,
    'networks': paths.TransitionNetwork,
    'cvs': paths.CollectiveVariable
}


class OPSStorage(storage.GeneralStorage):
    def __init__(self, backend, schema, class_info, fallbacks=None):
        # TODO: this will change to match the current notation
        super(OPSStorage, self).__init__(backend, schema, class_info,
                                         fallbacks)
        self.n_snapshot_types = 0

    @classmethod
    def from_backend(cls, backend, schema=None, class_info=None,
                     simulation_classes=None, fallbacks=None):
        obj = cls.__new__(cls)
        schema = tools.none_to_default(schema, ops_schema)
        class_info = tools.none_to_default(class_info, ops_class_info)
        simulation_classes = tools.none_to_default(simulation_classes,
                                                   ops_simulation_classes)
        super(OPSStorage, obj).__init__(
            backend=backend,
            schema=schema,
            class_info=class_info,
            simulation_classes=simulation_classes,
            fallbacks=fallbacks
        )
        obj.n_snapshot_types = 0
        return obj

    def register_from_tables(self, table_names, classes):
        lookups = {}
        table_to_class = {tbl: cls for tbl, cls in zip(table_names, classes)}
        for table in table_names:
            logger.info("Attempting to register missing table {} ({})"\
                        .format(table, str(table_to_class[table])))
            if issubclass(table_to_class[table], paths.BaseSnapshot):
                lookups.update(snapshots.snapshot_registration_from_db(
                    storage=self,
                    schema=self.schema,
                    class_info=self.class_info,
                    table_name=table
                ))
        logger.info("Found {} possible lookups".format(len(lookups)))
        logger.info("Lookups for tables: " + str(lookups.keys()))
        class_info_list = [ClassInfo(table=table,
                                     cls=table_to_class[table],
                                     lookup_result=lookups[table])
                           for table in lookups]
        for info in class_info_list:
            info.set_defaults(self.schema)
            self.class_info.add_class_info(info)


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

