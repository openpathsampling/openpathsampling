"""

@author: JD Chodera
@author: JH Prinz
"""

from openpathsampling.engines import BaseSnapshot
import openpathsampling.engines.features as feats


@feats.attach_features([
    feats.state,
    feats.engine
])
class MSMSnapshot(BaseSnapshot):
    """
    Simulation snapshot. Only references to coordinates and velocities
    """
