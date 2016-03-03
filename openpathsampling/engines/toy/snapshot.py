"""

@author: JD Chodera
@author: JH Prinz
"""

from openpathsampling.engines import BaseSnapshot
import features


@features.base.attach_features([
    features.velocities,
    features.coordinates,
    features.xyz,
    features.topology
])
class ToySnapshot(BaseSnapshot):
    """
    Simulation snapshot. Only references to coordinates and velocities
    """
