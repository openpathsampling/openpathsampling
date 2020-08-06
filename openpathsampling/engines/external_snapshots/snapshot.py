import time

from openpathsampling.engines import features
from openpathsampling.engines.snapshot import BaseSnapshot
from . import features as ext_features

import logging
logger = logging.getLogger(__name__)


@features.base.attach_features([
    features.engine,
    features.coordinates,
    features.velocities,
    features.box_vectors,
    ext_features.file_info
])
class InternalizedMDSnapshot(BaseSnapshot):
    """
    Internalized version of standard external MD snapshot.

    This can be used, for example, to store intial conditions in an OPS
    storage file.
    """
    pass

@features.base.attach_features([
    features.engine,
    ext_features.coordinates,
    ext_features.velocities,
    ext_features.box_vectors,
    ext_features.file_info
])
class ExternalMDSnapshot(BaseSnapshot):
    """
    Snapshot for external MD engines

    Internally, this only stores the file_name and the file_position. All
    specific details (positions, velocities, box vectors) are loaded from
    file when requested.

    Parameters
    ----------
    file_name : str
        the name of the external file where the positions/velocities/etc.
        reside
    file_position : int
        position within the file; the engine should be able to load data for
        this specific snapshot based on this number
    engine : :class:`.DynamicsEngine`
        the engine associated with this snapshot
    """
    def __init__(self, file_name=None, file_position=None, engine=None):
        # these are done in place of calling super
        self._reversed = None
        self.__uuid__ = self.get_uuid()
        # these are the requried attributes
        self.file_name = file_name
        self.file_position = file_position
        self.engine = engine
        self.velocity_direction = 1  # by default; reversed flips it
        # these are containers for temporary data
        self._xyz = None
        self._velocities = None
        self._box_vectors = None

        self._internalized = None

    def load_details(self):
        """Cache coords, velocities, box vectors from the external file"""
        try:
            (xyz, vel, box) = self.engine.read_frame_data(
                self.file_name,
                self.file_position
            )
        except IndexError:
            # Out of bounds on buffer access (axis 0)
            logger.debug("Exception reading from %s[%d]", self.file_name,
                         self.file_position)
            time.sleep(self.engine.sleep_ms/10000.0)  # 1/10 the normal
            self.load_details()
        except RecursionError:
            raise RuntimeError("Unrecoverable error in load_details")
        else:
            self._xyz = xyz
            self._velocities = vel
            self._box_vectors = box

    def set_details(self, xyz, velocities, box_vectors):
        """Set coords, velocities, and box vectors.

        This is mainly used if OPS must modify/create a snapshot.

        Parameters
        ----------
        xyz : np.array
            unitless coordinates
        velocities : np.array
            velocities
        box_vectors : np.array
            unit cell for the periodic box
        """
        try:
            self.load_details()
        except:
            pass
        else:
            raise RuntimeError("Can't set details if frame already exists.")
        finally:
            self._xyz = xyz
            self._velocities = velocities
            self._box_vectors = box_vectors

    def clear_cache(self):
        """Remove internal details from snapshot.

        These details should always be accessible later using
        :meth:`.load_details`. Removing them allows them memory to be
        freed.
        """
        self._xyz = None
        self._velocities = None
        self._box_vectors = None

    def __repr__(self):
        num_str = "file_name=" + str(self.file_name)
        pos_str = "file_position=" + str(self.file_position)
        eng_str = "engine=" + repr(self.engine)
        args = ", ".join([num_str, pos_str, eng_str])
        return "{cls_str}(".format(cls_str=self.cls) + args + ")"

    def internalize(self):
        """Return a version of this snapshot with storable details.

        This allows these snapshots to be stored internally in OPS storage
        files, instead of only in external files. This is convenient to
        avoid the need to transfer files to remote computers.
        """
        if self._internalized is None:
            self._internalized = self.engine.InternalizedSnapshotClass(
                coordinates=self.coordinates,
                velocities=self.velocities,
                box_vectors=self.box_vectors,
                file_name=self.file_name,
                file_position=self.file_position,
                engine=self.engine.internalized_engine
            )
        return self._internalized
