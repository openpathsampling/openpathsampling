import numpy as np
from shared import ConfigurationStore

attributes = ['configuration', 'box_vectors', 'md', 'coordinates', 'xyz']
lazy = ['configuration']


def netcdfplus_init(store):
    store.storage.create_store('configurations', ConfigurationStore())

    store.create_variable(
        'configuration',
        'lazyobj.configurations',
        description="the snapshot index (0..n_configuration-1) of snapshot '{idx}'.",
        chunksizes=(1,)
    )


def coordinates(snapshot):
    """
    Returns
    -------
    coordinates: numpy.ndarray, shape=(atoms, 3), dtype=numpy.float32
        the atomic coordinates of the configuration. The coordinates are wrapped in a
        simtk.unit.Unit.
    """

    if snapshot.configuration is not None:
        return snapshot.configuration.coordinates

    return None


def box_vectors(snapshot):
    """
    Returns
    -------
    box_vectors: numpy.ndarray, shape=(3, 3), dtype=numpy.float32
        the box_vectors of the configuration. The coordinates are wrapped in a
        simtk.unit.Unit.
    """
    if snapshot.configuration is not None:
        return snapshot.configuration.box_vectors

    return None


def md(snapshot):
    """
    Returns
    -------
    md : mdtraj.Trajectory
        the actual trajectory object. Can be used with all functions from mdtraj

    Notes
    -----
    Rather slow since the topology has to be made each time. Try to avoid it
    """

    if snapshot.configuration is not None:
        n_atoms = snapshot.coordinates.shape[0]

        output = np.zeros([1, n_atoms, 3], np.float32)
        output[0, :, :] = snapshot.coordinates

        return md.Trajectory(output, snapshot.topology.md)
