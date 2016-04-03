"""

@author: JD Chodera
@author: JH Prinz
"""

from openpathsampling.engines import BaseSnapshot, SnapshotFactory
import features


# @feats.attach_features([
#     feats.velocities,
#     feats.coordinates,
#     feats.xyz,
#     feats.topology
# ])
# class ToySnapshot(BaseSnapshot):
#     """
#     Simulation snapshot. Only references to coordinates and velocities
#     """


# The following code does the same as above

ToySnapshot = SnapshotFactory(
    name='ToySnapshot',
    features=[
        features.velocities,
        features.coordinates,
        features.topology
    ],
    description="Simulation snapshot. Only references to coordinates and velocities",
    base_class=BaseSnapshot
)