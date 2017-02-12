from openpathsampling.engines import Topology


class ToyTopology(Topology):
    """
    Parameters
    ----------
    n_spatial : int
        the number of spatial degrees of freedom
    masses : numpy.ndarray (n_atoms, dtype=float)
        The masses associated with each atom
    pes : :class:`.PES`
        potential energy surface for this system
    n_atoms : int
        number of atoms (default is 1)
    """

    def __init__(self, n_spatial, masses, pes, n_atoms=1):
        super(ToyTopology, self).__init__(n_atoms=n_atoms, n_spatial=n_spatial)
        self.masses = masses
        self.pes = pes
