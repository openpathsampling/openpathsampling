from openpathsampling.netcdfplus import ObjectStore
from openpathsampling.tools import simtk_units_from_md_snapshot

from openpathsampling.snapshot_content import Momentum

import numpy as np

attributes = ['momentum', 'velocities', 'is_reversed']
lazy = ['momentum']
flip = ['is_reversed']

def netcdfplus_init(store):
    store.storage.create_store('momenta', MomentumStore())

    store.create_variable('momentum', 'lazyobj.momenta',
                        description="the snapshot index (0..n_momentum-1) 'frame' of snapshot '{idx}'.",
                        )

    store.create_variable('is_reversed', 'bool',
                        description="the indicator if momenta should be flipped.",
                        )


def velocities(self):
    """
    The velocities in the configuration. If the snapshot is reversed a
    copy of the original (unreversed) velocities is made which is then
    returned
    """
    if self.momentum is not None:
        if self.is_reversed:
            return -1.0 * self.momentum.velocities
        else:
            return self.momentum.velocities

    return None


class MomentumStore(ObjectStore):
    """
    An ObjectStore for Momenta. Allows to store Momentum() instances in a netcdf file.
    """

    def __init__(self):
        super(MomentumStore, self).__init__(Momentum, json=False)

    def to_dict(self):
        return {}

    def _save(self, momentum, idx):
        self.vars['velocities'][idx, :, :] = momentum.velocities

    def _load(self, idx):
        velocities = self.vars['velocities'][idx]

        momentum = Momentum(velocities=velocities)
        return momentum

    def velocities_as_numpy(self, frame_indices=None, atom_indices=None):
        """
        Return a block of stored velocities in the database as a numpy array.

        Parameters
        ----------
        frame_indices : list of int or None
            the indices of Momentum objects to be retrieved from the database.
            If `None` is specified then all indices are returned!
        atom_indices : list of int of None
            if not None only the specified atom_indices are returned. Might
            speed up reading a lot.
        """

        if frame_indices is None:
            frame_indices = slice(None)

        if atom_indices is None:
            atom_indices = slice(None)

        return self.variables['velocities'][frame_indices, atom_indices, :].astype(np.float32).copy()

    def velocities_as_array(self, frame_indices=None, atom_indices=None):
        """
        Returns a numpy array consisting of all velocities at the given indices

        Parameters
        ----------
        frame_indices : list of int
            momenta indices to be loaded
        atom_indices : list of int
            selects only the atoms to be returned. If None (Default) all atoms
            will be selected


        Returns
        -------
        numpy.ndarray, shape = (l,n)
            returns an array with `l` the number of frames and `n` the number
            of atoms
        """

        return self.velocities_as_numpy(frame_indices, atom_indices)

    def _init(self):
        """
        Initializes the associated storage to index momentums in it
        """

        super(MomentumStore, self)._init()

        n_atoms = self.storage.n_atoms
        n_spatial = self.storage.n_spatial

        units = simtk_units_from_md_snapshot(self.storage._template)

        self.create_variable('velocities', 'numpy.float32',
                           dimensions=('atom', 'spatial'),
                           description="the velocity of atom 'atom' in dimension " +
                                       "'coordinate' of momentum 'momentum'.",
                           chunksizes=(1, n_atoms, n_spatial),
                           simtk_unit=units['velocity']
                           )