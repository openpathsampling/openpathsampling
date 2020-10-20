from ..simstore import storage
from ..simstore import sql_backend

from ..simstore.serialization_helpers import to_json_obj as json_serializer
from ..simstore.serialization_helpers import from_json_obj as deserialize_sim
from ..simstore.serialization_helpers import import_class
from ..simstore.serialization_helpers import get_uuid, set_uuid
from ..simstore.serialization_helpers import default_find_uuids

from ..simstore.class_lookup import ClassIsSomething

from ..simstore import (
    StorableFunction, StorableFunctionResults, storable_function_find_uuids
)

import openpathsampling as paths
from openpathsampling.netcdfplus import StorableObject

from ..simstore import tools

from ..simstore.custom_json import (
    JSONSerializerDeserializer,
    numpy_codec, bytes_codec, uuid_object_codec,
)

from ..simstore import CallableCodec

from ..simstore.serialization import (
    ToDictSerializer, SchemaSerializer, SchemaDeserializer
)

from ..simstore.class_info import ClassInfo, ClassInfoContainer

from ..simstore import SQLStorageBackend  # TODO: generalize

from . import snapshots
import logging
logger = logging.getLogger(__name__)

# this defines the schema for data objects
ops_schema = {
    'samples': [('trajectory', 'uuid'), #'lazy'),  # TODO: JHP's CVs fail
                ('ensemble', 'uuid'),
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
    'steps': [('change', 'uuid'), ('active', 'uuid'), #('previous', 'lazy'),
              ('simulation', 'uuid'), ('mccycle', 'int')],
    'details': [('json', 'json_obj')],
    'storable_functions': [('json', 'json_obj')],#, ('class_idx', 'int')],
    'simulation_objects': [('json', 'json_obj')],#, ('class_idx', 'int')]
    'storage_objects': [('json', 'json_obj')]
}

# this includes any sql-specific metadata
ops_schema_sql_metadata = {}

# this defines the simulation object serializer for OPS
CODECS = [numpy_codec, bytes_codec, uuid_object_codec]

UNSAFE_CODECS = CODECS + [CallableCodec()]
SAFE_CODECS = CODECS + [CallableCodec({'safemode': True})]

class MoveChangeDeserializer(SchemaDeserializer):
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



class OPSSpecialLookup(object):
    """Separate object to handle special lookups

    This is separate because, in addition to the functionality it encodes,
    it also acts as a cache to reduce the number of isinstance calls (and
    hopfully speed up the identification of types)
    """
    # TODO: overall lookup should be like this: all lookups map a class to a
    # lookup function -- often that lookup function just returns the class
    # TODO: especially annoying to have StorableFunction as a special
    special_superclasses = (paths.BaseSnapshot, paths.MoveChange,
                            paths.Details, StorableFunction,
                            storage.GeneralStorage, SQLStorageBackend)
    snapshot_lookup_function = \
            lambda self, snap: (get_uuid(snap.engine), snap.__class__)
    details_lookup_function = lambda self, details: paths.Details
    movechange_lookup_function = lambda self, change: paths.MoveChange
    sf_lookup_function = lambda self, func: StorableFunction
    storage_obj_lookup_function = lambda self, func: storage.GeneralStorage

    def __init__(self):
        self.secondary_lookups = {}
        # self.special_classes = set()
        # self.non_special_classes = set()
        is_special_func = lambda obj: \
                isinstance(obj, self.special_superclasses)
        self.is_special = ClassIsSomething(is_special_func)

    def __call__(self, item):
        cls = item.__class__
        if cls in self.secondary_lookups:
            return self.secondary_lookups[cls](item)

        if isinstance(item, paths.BaseSnapshot):
            self.secondary_lookups[cls] = self.snapshot_lookup_function
        elif isinstance(item, StorableFunction):
            self.secondary_lookups[cls] = self.sf_lookup_function
        elif isinstance(item, (storage.GeneralStorage, SQLStorageBackend)):
            self.secondary_lookups[cls] = self.storage_obj_lookup_function
        elif isinstance(item, paths.MoveChange):
            self.secondary_lookups[cls] = self.movechange_lookup_function
        elif isinstance(item, paths.Details):
            # TODO: this should be removed, since all Details classes are
            # equivalent -- unfortunately, JHP's LoaderProxy breaks in that
            # case
            self.secondary_lookups[cls] = self.details_lookup_function

        return self.secondary_lookups[cls](item)

class OPSClassInfoContainer(ClassInfoContainer):
    def __init__(self, default_info, sfr_info=None, schema=None,
                 class_info_list=None):
        super(OPSClassInfoContainer, self).__init__(default_info,
                                                    sfr_info,
                                                    schema,
                                                    class_info_list)
        self.n_snapshot_types = 0
        self.special_lookup_object = OPSSpecialLookup()

    def is_special(self, item):
        return self.special_lookup_object.is_special(item)

    def special_lookup_key(self, item):
        return self.special_lookup_object(item)

    def add_missing_table_from_instance(self, lookup, obj):
        if isinstance(obj, paths.BaseSnapshot):
            schema, class_info_list = snapshots.snapshot_registration_info(
                obj, self.n_snapshot_types
            )
            schema = snapshots.replace_schema_dimensions(schema,
                                                         obj.engine.descriptor)
            self.register_info(class_info_list, schema)
            self.n_snapshot_types += 1

unsafe_ops_codecs = JSONSerializerDeserializer(UNSAFE_CODECS)
safe_ops_codecs = JSONSerializerDeserializer(SAFE_CODECS)

def _build_ops_serializer(schema, safe_codecs, unsafe_codecs):
    # TODO: why is this using deserialize_sim instead of the codec
    # deserializer?  probably need to change that for safemode
    ops_class_info = OPSClassInfoContainer(
        default_info=ClassInfo(
            table='simulation_objects',
            cls=StorableObject,
            serializer=unsafe_codecs.simobj_serializer,
            deserializer=unsafe_codecs.simobj_deserializer,
            safe_deserializer=safe_codecs.simobj_serializer,
            find_uuids=default_find_uuids
        ),
        sfr_info=ClassInfo(
            table="function_results",
            cls=StorableFunctionResults,
            serializer=StorableFunctionResults.to_dict,
            deserializer=lambda x: x,  # deserializer not used
            find_uuids=default_find_uuids
        ),
        schema=schema,
        class_info_list=[
            ClassInfo(table='samples', cls=paths.Sample),
            ClassInfo(table='sample_sets', cls=paths.SampleSet),
            ClassInfo(table='trajectories', cls=paths.Trajectory),
            ClassInfo(table='move_changes', cls=paths.MoveChange,
                      deserializer=MoveChangeDeserializer(
                          schema=schema,
                          table='move_changes'
                      )),
            ClassInfo(table='steps', cls=paths.MCStep),
            ClassInfo(table='details', cls=paths.Details,
                      serializer=safe_codecs.simobj_serializer,
                      deserializer=safe_codecs.simobj_deserializer),
            ClassInfo(table='storable_functions',
                      cls=StorableFunction,
                      find_uuids=storable_function_find_uuids,
                      serializer=unsafe_codecs.simobj_serializer,
                      deserializer=unsafe_codecs.simobj_deserializer),
            ClassInfo(table="storage_objects",
                      cls=storage.GeneralStorage,
                      serializer=safe_codecs.simobj_serializer,
                      deserializer=safe_codecs.simobj_deserializer,
                      find_uuids=default_find_uuids),
        ]
    )

    for info in ops_class_info.class_info_list:
        info.set_defaults(schema)

    return ops_class_info

ops_class_info = _build_ops_serializer(ops_schema, safe_ops_codecs,
                                       unsafe_ops_codecs)

# this will create the pseudo-tables used to find specific objects
ops_simulation_classes = {
    'volumes': paths.Volume,
    'ensembles': paths.Ensemble,
    'pathsimulators': paths.PathSimulator,
    'pathmovers': paths.PathMover,
    'networks': paths.TransitionNetwork,
    'cvs': (paths.CollectiveVariable, StorableFunction),
    'engines': paths.engines.DynamicsEngine
}  # TODO: add more to these


class Storage(storage.GeneralStorage):
    def __init__(self, backend, schema, class_info, fallbacks=None,
                 safemode=False):
        # TODO: this will change to match the current notation
        super(Storage, self).__init__(backend, schema, class_info,
                                      fallbacks, safemode)

        self.n_snapshot_types = 0

    def sync_all(self):
        self.save(self._stashed)
        self._stashed = []

    @classmethod
    def from_backend(cls, backend, schema=None, class_info=None,
                     simulation_classes=None, fallbacks=None,
                     safemode=False):
        # quick exit if this storage is known
        exists = None
        if backend.identifier[1] != 'w':
            exists = cls._known_storages.get(backend.identifier, None)
        if exists is not None:
            return exists
        obj = cls.__new__(cls)
        schema = tools.none_to_default(schema, ops_schema)
        class_info = tools.none_to_default(class_info, ops_class_info)
        simulation_classes = tools.none_to_default(simulation_classes,
                                                   ops_simulation_classes)
        super(Storage, obj).__init__(
            backend=backend,
            schema=schema,
            class_info=class_info,
            simulation_classes=simulation_classes,
            fallbacks=fallbacks
        )
        obj.n_snapshot_types = 0
        return obj

    def to_dict(self):
        return {'backend': self.backend,
                'fallbacks': self.fallbacks,
                'safemode': self.safemode}

    @classmethod
    def from_dict(cls, dct):
        return cls.from_backend(**dct)

    def __reduce__(self):
        return (self.from_dict, (self.to_dict(),))

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
        self.class_info.register_info(class_info_list, self.schema)
        # for info in class_info_list:
            # info.set_defaults(self.schema)
            # self.class_info.add_class_info(info)


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

