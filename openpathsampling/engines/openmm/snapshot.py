"""

@author: JD Chodera
@author: JH Prinz
"""

from openpathsampling.engines import BaseSnapshot, SnapshotFactory
import features


@features.base.attach_features([
    features.velocities,
    features.coordinates,
    features.box_vectors,
    features.engine
])
class MDSnapshot(BaseSnapshot):
    """
    A fast MDSnapshot
    """

# The following code does the same as above
#
#  MDSnapshot = SnapshotFactory(
#     name='MDSnapshot',
#     features=[
#         features.velocities,
#         features.coordinates,
#         features.box_vectors,
#         features.engine
#     ],
#     description="A fast MDSnapshot",
#     base_class=BaseSnapshot
# )


@features.base.attach_features([
    features.statics,
    features.kinetics,
    features.engine
])
class Snapshot(BaseSnapshot):
    """
    The standard MDSnapshot supporting coordinate, velocities and box_vectors
    """

    StaticContainer = features.StaticContainer
    KineticContainer = features.KineticContainer

    @staticmethod
    def construct(coordinates=None, box_vectors=None, velocities=None, engine=None):
        """
        Construct a new snapshot from numpy arrays

        This will create the container objects and return a Snapshot object. Mostly a helper
        to allow for easier creation.

        Parameters
        ----------
        coordinates : numpy.array, shape = (atoms, spatial)
            the atomic coordinates
        box_vectors : numpy.array, shape = (spatial, spatial)
            the box vectors
        velocities : numpy.array, shape = (atoms, spatial)
            the atomic velocities

        Returns
        -------
        :obj:`Snapshot`
            the created `Snapshot` object
        """
        statics = Snapshot.StaticContainer(coordinates=coordinates, box_vectors=box_vectors)
        kinetics = Snapshot.KineticContainer(velocities=velocities)

        return Snapshot(engine=engine, statics=statics, kinetics=kinetics)

    @property
    def topology(self):
        return self.engine.topology