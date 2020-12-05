"""

@author: JD Chodera
@author: JH Prinz
"""

from openpathsampling.engines import BaseSnapshot, SnapshotFactory
from . import features


@features.base.attach_features([
    features.velocities,
    features.coordinates,
    features.box_vectors,
    features.engine
])
class MDSnapshot(BaseSnapshot):
    """
    A fast MD snapshot, which does not proxy the coordinates/velocities.
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
    features.masses,
    features.instantaneous_temperature,
    features.engine,
    features.traj_quantities,
])
class Snapshot(BaseSnapshot):
    """
    The standard snapshot for MD, based on statics and kinetics proxies.
    """

    StaticContainer = features.StaticContainer
    KineticContainer = features.KineticContainer

    @staticmethod
    def construct(
            coordinates=None,
            box_vectors=None,
            velocities=None,
            statics=None,
            kinetics=None,
            engine=None):
        """
        Construct a new snapshot from numpy arrays

        This will create the container objects and return a Snapshot object.
        Mostly a helper to allow for easier creation.

        You can either use coordinates and velocities and/or statics and
        kinetics objects. If both are present the more complex (statics
        and kinetics) will be used

        Parameters
        ----------
        coordinates : numpy.array, shape ``(n_atoms, n_spatial)``
            the atomic coordinates
        box_vectors : numpy.array, shape ``(n_spatial, n_spatial)``
            the box vectors
        velocities : numpy.array, shape ``(n_atoms, n_spatial)``
            the atomic velocities
        statics : `openpathsampling.engines.openmm.StaticContainer`
            the statics container if it already exists
        kinetics : `openpathsampling.engines.openmm.KineticContainer`
            the kinetics container if it already exists

        engine : :obj:`openpathsampling.engines.DynamicsEngine`
            the engine that should be referenced as the one used to
            generate the object

        Returns
        -------
        :obj:`Snapshot`
            the created `Snapshot` object
        """
        if statics is None:
            statics = Snapshot.StaticContainer(
                coordinates=coordinates,
                box_vectors=box_vectors,
                engine=engine
            )

        if kinetics is None:
            kinetics = Snapshot.KineticContainer(
                velocities=velocities,
                engine=engine
            )

        return Snapshot(
            engine=engine,
            statics=statics,
            kinetics=kinetics
        )

    @property
    def topology(self):
        return self.engine.topology
