from openpathsampling.engines import Topology


class HOOMDTopology(Topology):
    """
    Parameters
    ----------
    n_spatial : int
        Number of spatial degrees of freedom.
    n_atoms : int
        Number of particles.
    """

    def __init__(self, n_spatial, n_atoms=1):
        super(HOOMDTopology, self).__init__(n_atoms=n_atoms, n_spatial=n_spatial)
