from collections import abc
import itertools

from ..simstore import StorageTable

# we need a separate class for the snapshots table because the snapshots
# table actually combines multiple tables
class SnapshotsTable(abc.Sequence):
    def __init__(self, storage):
        self.storage = storage
        self.tables = {}
        self.update_tables()

    def update_tables(self):
        backend = self.storage.backend
        snapshot_table_names = [table for table in backend.table_to_class
                                if table.startswith('snapshot')]
        new_tables = [table for table in snapshot_table_names
                      if table not in self.tables]

        self.tables.update({table: StorageTable(self.storage, table)
                            for table in new_tables})

    def __iter__(self):
        return itertools.chain.from_iterable(self.tables.values())

    def __getitem__(self, item):
        if len(self.tables) == 0:
            raise IndexError("index out of range")

        if item < 0:
            item += len(self)

        count = 0
        for table in self.tables.values():
            len_table = len(table)
            count += len_table
            if count > item:
                count -= len_table
                break

        return table[item - count]

    def __len__(self):
        return sum(len(table) for table in self.tables.values())

    def snapshots_for_engine(self, engine):
        engine_dict = {table[0].engine: table
                       for table in self.tables.values()}
        return engine_dict[engine]

    def cache_all(self):
        # unit testing this will be hard; prob requires real storage
        old_blocksize = self.iter_block_size
        self.iter_block_size = len(self)
        _ = list(iter(self))
        self.iter_block_size = old_blocksize

    def save(self, obj):
        self.storage.save(obj)
