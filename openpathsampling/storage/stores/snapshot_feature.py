import logging

from .snapshot_base import BaseSnapshotStore

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


class FeatureSnapshotStore(BaseSnapshotStore):
    """
    An ObjectStore for Snapshots in netCDF files.
    """

    def __init__(self, descriptor):
        super(FeatureSnapshotStore, self).__init__(descriptor)

    @property
    def classes(self):
        return self.snapshot_class.__features__.classes

    @property
    def storables(self):
        return self.snapshot_class.__features__.storables

    def _set(self, idx, snapshot):
        [self.write(attr, idx, snapshot) for attr in self.storables]

    def _get(self, idx, snapshot):
        [setattr(snapshot, attr, self.vars[attr][idx])
         for attr in self.storables]

    def initialize(self):
        super(FeatureSnapshotStore, self).initialize()

        for dim, size in self._dimensions.items():
            self.storage.create_dimension(self.prefix + dim, size)

        for feature in self.classes:
            if hasattr(feature, 'netcdfplus_init'):
                feature.netcdfplus_init(self)

        self.storage.sync()
