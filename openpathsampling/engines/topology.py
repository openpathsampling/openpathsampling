from openpathsampling.netcdfplus import StorableNamedObject


class Topology(StorableNamedObject):
    """
    Topology is the object that contains all information about the structure
    of the system to be simulated.

    Attributes
    ----------
    n_atoms : int
        number of atoms
    n_spatial : int
        number of spatial dimensions, default is 3
    """

    def __init__(self, n_atoms, n_spatial=3):
        super(Topology, self).__init__()
        self.n_atoms = n_atoms
        self.n_spatial = n_spatial
