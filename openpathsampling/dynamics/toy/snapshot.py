"""

@author: JD Chodera
@author: JH Prinz
"""

from openpathsampling.dynamics import FeatureSnapshot
import features


@features.base.set_features(
    features.velocities,
    features.coordinates,
    features.xyz,
    features.topology
)
class ToySnapshot(FeatureSnapshot):
    """
    Simulation snapshot. Only references to coordinates and velocities
    """
