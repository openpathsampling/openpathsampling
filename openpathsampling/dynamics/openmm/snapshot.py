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