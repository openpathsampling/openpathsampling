"""

@author: JD Chodera
@author: JH Prinz
"""

from openpathsampling.dynamics import BaseSnapshot
import features


@features.base.attach_features([
    features.velocities,
    features.coordinates,
    features.box_vectors,
    features.xyz,
    features.topology
])
class MDSnapshot(BaseSnapshot):
    """
    A fast MDSnapshot
    """


@features.base.attach_features([
    features.configuration,
    features.momentum,
    features.xyz,
    features.topology  # for compatibility
])
class Snapshot(BaseSnapshot):
    """
    The standard MDSnapshot supporting coordinate, velocities and box_vectors
    """

    Configuration = features.Configuration
    Momentum = features.Momentum

    @staticmethod
    def construct(coordinates=None, box_vectors=None, velocities=None, topology=None):
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
        configuration = Snapshot.Configuration(coordinates=coordinates, box_vectors=box_vectors)
        momentum = Snapshot.Momentum(velocities=velocities)

        return Snapshot(topology=topology, configuration=configuration, momentum=momentum)