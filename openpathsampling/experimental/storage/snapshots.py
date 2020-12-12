# NOTE: this is part of the OPS-specific stuff
from ..simstore.class_info import ClassInfo
from ..simstore.uuids import get_uuid
from ..simstore.proxy import GenericLazyLoader

def _nested_schema_entries(schema_entries, lazies):
    """Recursive algorithm to create all schema entries
    """
    entries = []
    schema = {}
    for (feat_name, feat_type) in schema_entries:
        if not isinstance(feat_type, str):
            loc_entr, loc_schema = _nested_schema_entries(feat_type, lazies)
            schema.update(loc_schema)
            schema.update({feat_name: loc_entr})
            feat_type = 'lazy' if feat_name in lazies else 'uuid'
        entries.append((feat_name, feat_type))
    return entries, schema


def schema_from_entries(features, lazies):
    """Build the schema dict from the features.

    Note that the resulting dict has two types of placeholders, compared to
    the actual snapshot: the snapshot name will be changed by the storage,
    and and dimensions in the type definitions will be replaced by values
    from the SnapshotDescriptor.

    Parameters
    ----------
    features : list of modules
        the feature modules to be used
    lazies : list of str
        the feature names that have been marked as lazy

    Returns
    -------
    dict
        the schema dictionary, ready for replacement of dimensions by
        SnapshotDescriptor
    """
    # load entries from all features, recurse over them to find all
    # subentries, and then return the dict tht comes from it all
    schema_entries = sum([feat.schema_entries for feat in features
                          if hasattr(feat, 'schema_entries')], [])
    entries, schema = _nested_schema_entries(schema_entries, lazies)
    schema.update({'snapshot': entries})
    return schema


def schema_for_snapshot(snapshot):
    return schema_from_entries(features=snapshot.__features__.classes,
                               lazies=snapshot.__features__.lazy)


def replace_schema_dimensions(schema, descriptor):
    descriptor_dict = {desc[0]: desc[1] for desc in descriptor
                       if desc[0] != 'class'}
    for (table, entries) in schema.items():
        schema[table] = [
            (attr, type_name.format(**descriptor_dict))
            for (attr, type_name) in entries
        ]
    return schema


def snapshot_registration_from_db(storage, schema, class_info, table_name):
    # TODO: snapshot tables always have `snapshotNUM`; this should be used
    # to identify other related tables in the DB
    cls = storage.backend.table_to_class[table_name]
    representative_row = storage.backend.get_representative(table_name)
    engine_uuid = representative_row.engine
    lookup_result = (engine_uuid, cls)
    proposed_lookups = {table_name: lookup_result}
    attributes = schema[table_name]
    # for (attr, type_name) in attributes:
        # is_object = type_name in ['lazy', 'uuid', 'uuid_list']
        # is_table = attr in schema
        # if is_object and is_table:
            # cls = storage.backend.table_to_class[attr]
            # proposed_lookups[attr] = (engine_uuid, cls)
    return proposed_lookups


def snapshot_registration_info(snapshot_instance, snapshot_number):
    schema = schema_for_snapshot(snapshot_instance)
    real_table = {table: table + str(snapshot_number) for table in schema}
    real_schema = {real_table[table]: entries
                   for (table, entries) in schema.items()}
    engine = snapshot_instance.engine
    engine_uuid = get_uuid(engine)
    snapshot_info = ClassInfo(table=real_table['snapshot'],
                              cls=snapshot_instance.__class__,
                              lookup_result=(engine_uuid,
                                             snapshot_instance.__class__))
    attr_infos = []
    for table in [tbl for tbl in schema.keys() if tbl != 'snapshot']:
        obj = getattr(snapshot_instance, table)
        if isinstance(obj, GenericLazyLoader):
            obj = obj.load()
        attr_infos.append(ClassInfo(table=real_table[table],
                                    cls=obj.__class__,
                                    lookup_result=(engine_uuid,
                                                   obj.__class__)))
    class_info_list = [snapshot_info] + attr_infos
    return real_schema, class_info_list
