import collections

from ..simstore import storage
from ..simstore import sql_backend

from ..simstore.attribute_handlers import DEFAULT_HANDLERS

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
    JSONSerializerDeserializer, DEFAULT_CODECS
)

from ..simstore import CallableCodec

from ..simstore.serialization import (
    ToDictSerializer, SchemaSerializer, SchemaDeserializer
)

from ..simstore.class_info import ClassInfo, ClassInfoContainer

from ..simstore import SQLStorageBackend  # TODO: generalize

from . import snapshots
from .snapshots_table import SnapshotsTable
from .collective_variables import CollectiveVariable

try:
    from .simtk_unit import simtk_quantity_codec, SimtkQuantityHandler
except ImportError:
    simtk_quantity_codec = None
    SimtkQuantityHandler = None

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
EXTRA_CODECS = [simtk_quantity_codec] if simtk_quantity_codec else []
EXTRA_HANDLERS = [SimtkQuantityHandler] if SimtkQuantityHandler else []
CODECS = DEFAULT_CODECS + EXTRA_CODECS
HANDLERS = DEFAULT_HANDLERS + EXTRA_HANDLERS

UNSAFE_CODECS = CODECS + [CallableCodec()]
SAFE_CODECS = CODECS + [CallableCodec({'safemode': True})]

ENGINE_LOOKUP_CLASSES = (
    paths.BaseSnapshot,
    paths.engines.features.shared.KineticContainer,
    paths.engines.features.shared.StaticContainer,
)

class MoveChangeDeserializer(SchemaDeserializer):
    # in general, I think it would be better to reorg MoveChange to only be
    # one class, but this is aimed at fixing problems with reloading
    # MoveChange objects
    def __init__(self, schema, table, handlers):
        super(MoveChangeDeserializer, self).__init__(
            schema=schema,
            table=table,
            cls=None,
            handlers=[]
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


SpecialLookup = collections.namedtuple("SpecialLookup",
                                       ['cls', 'secondary_lookup'])

def engine_secondary_lookup(obj):
    return (get_uuid(obj.engine), obj.__class__)

SPECIAL_LOOKUPS = [
    SpecialLookup(paths.BaseSnapshot, engine_secondary_lookup),
    SpecialLookup(paths.engines.features.shared.KineticContainer,
                  engine_secondary_lookup),
    SpecialLookup(paths.engines.features.shared.StaticContainer,
                  engine_secondary_lookup),
    # most of the rest of these are just to map weird subclasses
    SpecialLookup(paths.MoveChange, lambda change: paths.MoveChange),
    SpecialLookup(paths.Details, lambda details: paths.Details),
    SpecialLookup(storage.GeneralStorage,
                  lambda st: storage.GeneralStorage),
    SpecialLookup(StorableFunction, lambda sf: StorableFunction),
]

class SpecialLookups(object):
    # TODO: looks like we can move Specials over to SimStore
    def __init__(self, lookups=None):
        lookups = tools.none_to_default(lookups, [])
        self.lookups = []
        self.superclass_secondary_lookups = {}
        self.special_superclasses = []
        self.secondary_lookups = {}
        self._is_special = None
        for lookup in lookups:
            self.register(lookup)

    def register(self, lookup):
        self.lookups.append(lookup)
        self.special_superclasses.append(lookup.cls)
        self.superclass_secondary_lookups[lookup.cls] = lookup.secondary_lookup
        self._is_special = ClassIsSomething(
            lambda obj: isinstance(obj, tuple(self.special_superclasses))
        )

    def is_special(self, obj):
        return self._is_special(obj)

    def __call__(self, obj):
        cls = obj.__class__
        if cls in self.secondary_lookups:
            return self.secondary_lookups[cls](obj)

        for lookup in self.lookups:
            if isinstance(obj, lookup.cls):
                self.secondary_lookups[cls] = lookup.secondary_lookup
                break

        return self.secondary_lookups[cls](obj)


class OPSClassInfoContainer(ClassInfoContainer):
    def __init__(self, default_info, sfr_info=None, schema=None,
                 class_info_list=None):
        super(OPSClassInfoContainer, self).__init__(default_info,
                                                    sfr_info,
                                                    schema,
                                                    class_info_list,
                                                    HANDLERS)
        self.n_snapshot_types = 0
        # self.special_lookup_object = OPSSpecialLookup()
        self.special_lookup_object = SpecialLookups(SPECIAL_LOOKUPS)

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
                          table='move_changes',
                          handlers=HANDLERS
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
        info.set_defaults(schema, HANDLERS)

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
    'schemes': paths.MoveScheme,
    'cvs': (paths.CollectiveVariable, CollectiveVariable),
    'engines': paths.engines.DynamicsEngine
}  # TODO: add more to these


class Storage(storage.GeneralStorage):
    def __init__(self, filename, mode='r', fallbacks=None, safemode=False):
        # TODO: this will change to match the current notation
        backend = sql_backend.SQLStorageBackend(filename, mode=mode)
        self.snapshots = None
        super(Storage, self).__init__(
            backend=backend,
            schema=ops_schema,
            class_info=ops_class_info,
            simulation_classes=ops_simulation_classes,
            fallbacks=fallbacks,
            safemode=safemode
        )

        self.snapshots = SnapshotsTable(self)
        self.snapshots.update_tables()

    @property
    def n_snapshot_types(self):
        return len(self.snapshots.tables)

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
        obj.snapshots = None
        super(Storage, obj).__init__(
            backend=backend,
            schema=schema,
            class_info=class_info,
            simulation_classes=simulation_classes,
            fallbacks=fallbacks
        )
        obj.snapshots = SnapshotsTable(obj)
        obj.snapshots.update_tables()
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

    @property
    def movechanges(self):
        return self.move_changes

    @property
    def samplesets(self):
        return self.sample_sets

    def register_from_tables(self, table_names, classes):
        lookups = {}
        table_to_class = {tbl: cls for tbl, cls in zip(table_names, classes)}
        for table in table_names:
            logger.info("Attempting to register missing table {} ({})"\
                        .format(table, str(table_to_class[table])))
            if issubclass(table_to_class[table], ENGINE_LOOKUP_CLASSES):
                lookups.update(snapshots.snapshot_registration_from_db(
                    storage=self,
                    schema=self.schema,
                    class_info=self.class_info,
                    table_name=table
                ))
                if self.snapshots is not None:
                    self.snapshots.update_tables()
        logger.info("Found {} possible lookups".format(len(lookups)))
        logger.info("Lookups for tables: " + str(lookups.keys()))
        # TODO: no, this is wrong; we need custom lookup
        class_info_list = [ClassInfo(table=table,
                                     cls=table_to_class[table],
                                     lookup_result=lookups[table])
                           for table in lookups]
        # self.class_info.register_info(class_info_list, self.schema)
        self.register_schema(self.schema, class_info_list, read_mode=True)
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
            self.snapshots.update_tables()  # increments snapshot types, too

