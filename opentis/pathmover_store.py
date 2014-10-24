from object_storage import ObjectStorage
from pathmover import PathMover, MoveDetails
from shooting import ShootingPoint, ShootingPointSelector

class PathMoverStorage(ObjectStorage):
    def __init__(self, storage):
        super(PathMoverStorage, self).__init__(storage, PathMover, named=True, json=True, identifier='json')

class ShootingPointStorage(ObjectStorage):
    def __init__(self, storage):
        super(ShootingPointStorage, self).__init__(storage, ShootingPoint, named=True, json=True)

class ShootingPointSelectorStorage(ObjectStorage):
    def __init__(self, storage):
        super(ShootingPointSelectorStorage, self).__init__(storage, ShootingPointSelector, named=True, json=True, identifier='json')

class MoveDetailsStorage(ObjectStorage):
    def __init__(self, storage):
        super(MoveDetailsStorage, self).__init__(storage, MoveDetails, named=False, json=True, identifier='json')