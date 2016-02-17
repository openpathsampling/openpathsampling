"""

@author: JD Chodera
@author: JH Prinz
"""

from openpathsampling.dynamics import FeatureSnapshot
import features


@features.base.set_features(
    features.velocities,
    features.coordinates,
    features.box_vectors,
    features.xyz,
    features.topology
)
class MDSnapshot(FeatureSnapshot):
    """
    A fast MDSnapshot
    """


@features.base.set_features(
    features.configuration,
    features.momentum,
    features.xyz,
    features.topology  # for compatibility
)
class Snapshot(FeatureSnapshot):
    """
    The standard MDSnapshot supporting coordinate, velocities and box_vectors
    """

    Configuration = features.Configuration
    Momentum = features.Momentum