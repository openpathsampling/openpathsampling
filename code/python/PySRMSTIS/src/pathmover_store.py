from object_storage import ObjectStorage, addcache
from pathmover import PathMover
from shooting import ShootingPoint, ShootingPointSelector

class PathMoverStorage(ObjectStorage):
    def __init__(self, storage):
        super(PathMoverStorage, self).__init__(storage, PathMover, named=True, json=True)

class ShootingPointStorage(ObjectStorage):
    def __init__(self, storage):
        super(ShootingPointStorage, self).__init__(storage, ShootingPoint, named=True, json=True)

class ShootingPointSelectorStorage(ObjectStorage):
    def __init__(self, storage):
        super(ShootingPointSelectorStorage, self).__init__(storage, ShootingPointSelector, named=True, json=True)