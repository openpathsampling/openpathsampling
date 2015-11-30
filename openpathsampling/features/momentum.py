from openpathsampling.netcdfplus import ObjectStore
from openpathsampling.tools import units_from_snapshot

from openpathsampling.snapshot import Momentum

import numpy as np

_variables = ['momentum']


def _init(store):
    store.storage.create_store('momenta', MomentumStore())

    store.create_variable('momentum', 'lazyobj.momenta',
                        description="the snapshot index (0..n_momentum-1) 'frame' of snapshot '{idx}'.",
                        chunksizes=(1,)
                        )


class MomentumStore(ObjectStore):
    """
    An ObjectStore for Momenta. Allows to store Momentum() instances in a netcdf file.
    """

    def __init__(self):
        super(MomentumStore, self).__init__(Momentum, json=False)

    def _save(self, momentum, idx):
        self.vars['velocities'][idx, :, :] = momentum.velocities

        if momentum.kinetic_energy is not None:
            self.vars['kinetic_energy'][idx] = momentum.kinetic_energy

    def _load(self, idx):
        velocities = self.vars['velocities'][idx]
        kinetic_energy = self.vars['kinetic_energy'][idx]

        momentum = Momentum(velocities=velocities, kinetic_energy=kinetic_energy)

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

        units = units_from_snapshot(self.storage._template)

        self.create_variable('velocities', 'numpy.float32',
                           dimensions=('atom', 'spatial'),
                           description="the velocity of atom 'atom' in dimension " +
                                       "'coordinate' of momentum 'momentum'.",
                           chunksizes=(1, n_atoms, n_spatial),
                           simtk_unit=units['velocity']
                           )

        self.create_variable('kinetic_energy', 'float',
                           chunksizes=(1,),
                           simtk_unit=units['energy']
                           )