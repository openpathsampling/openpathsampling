from .tools import group_by

def extract_backend_metadata(metadata, table, column):
    col_metadata = {}  # default
    if metadata:
        try:
            col_metadata = metadata[table][column]
        except KeyError:
            # doesn't matter whether missing key was table or column
            pass
    return col_metadata

def backend_table_list_consistency(backend):
    """Check that the stored list of tables matches the backend's schema"""
    schema_tables = backend.schema.keys()
    pass

def backend_table_inconsistencies(backend, table_name):
    """
    Check that UUIDs in a table match the UUIDs pointing to that table
    """
    all_uuids = list(backend.table_iterator('uuid'))
    all_uuids_by_table = group_by(all_uuids,
                                  lambda r: backend.uuid_row_to_table_name)
    tables = list(backend.table_iterator('tables'))
    pass
