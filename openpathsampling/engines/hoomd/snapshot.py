"""
@author: Bradley Dice
"""

from openpathsampling.engines import BaseSnapshot

from . import features


@features.base.attach_features(
    [
        features.statics,
        features.kinetics,
        features.engine,
    ]
)
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
        engine=None,
    ):
        """Construct a new snapshot from numpy arrays.

        This will create the container objects and return a Snapshot object.
        Mostly a helper to allow for easier creation.

        You can either use coordinates and velocities and/or statics and
        kinetics objects. If both are present, the more complex (statics
        and kinetics) will be used.

        Parameters
        ----------
        coordinates : numpy.array, shape ``(n_atoms, n_spatial)``
            The particle coordinates.
        box_vectors : numpy.array, shape ``(n_spatial, n_spatial)``
            The box vectors.
        velocities : numpy.array, shape ``(n_atoms, n_spatial)``
            The particle velocities.
        statics : `openpathsampling.engines.openmm.StaticContainer`
            The statics container if it already exists.
        kinetics : `openpathsampling.engines.openmm.KineticContainer`
            The kinetics container if it already exists.
        engine : `openpathsampling.engines.DynamicsEngine`
            The engine that should be referenced as the one used to
            generate the object.

        Returns
        -------
        `Snapshot`
            the created `Snapshot` object

        """
        if statics is None:
            statics = Snapshot.StaticContainer(
                coordinates=coordinates, box_vectors=box_vectors, engine=engine
            )

        if kinetics is None:
            kinetics = Snapshot.KineticContainer(velocities=velocities, engine=engine)

        return Snapshot(engine=engine, statics=statics, kinetics=kinetics)
