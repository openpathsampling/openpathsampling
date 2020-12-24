import numpy as np
from .shared import StaticContainerStore, StaticContainer, unmask_quantity
from openpathsampling.netcdfplus import WeakLRUCache
import openpathsampling as paths

variables = ['statics']
lazy = ['statics']

storables = ['statics']

dimensions = ['n_atoms', 'n_spatial']

_length_unit = "simtk(unit.nanometer)"
_array32 = "ndarray.float32"
schema_entries = [
    ('statics', [
        ('coordinates',
         '{length_unit}*{array32}({{n_atoms}},{{n_spatial}})'.format(
             length_unit=_length_unit, array32=_array32
        )),
        ('box_vectors',
         '{length_unit}*{array32}({{n_spatial}},{{n_spatial}})'.format(
             length_unit=_length_unit, array32=_array32
        )),
        ('engine', 'uuid'),
    ]),
]


def netcdfplus_init(store):
    static_store = StaticContainerStore()
    static_store.set_caching(WeakLRUCache(10000))

    name = store.prefix + 'statics'

    static_store.set_dimension_prefix_store(store)

    store.storage.create_store(name, static_store, False)

    store.create_variable(
        'statics',
        'lazyobj.' + name,
        description="the snapshot index (0..n_configuration-1) of "
                    "snapshot '{idx}'.")


@property
def coordinates(snapshot):
    """
    Returns
    -------
    coordinates: numpy.ndarray, shape=(atoms, 3), dtype=numpy.float32
        the atomic coordinates of the configuration. The coordinates are
        wrapped in a `simtk.unit.Unit`.
    """

    if snapshot.statics is not None:
        return unmask_quantity(snapshot.statics.coordinates)

    return None


@coordinates.setter
def coordinates(self, value):
    if value is not None:
        sc = StaticContainer(coordinates=value,
                             box_vectors=self.box_vectors,
                             engine=self.engine)
    else:
        sc = None

    self.statics = sc


@property
def box_vectors(snapshot):
    """
    Returns
    -------
    box_vectors: numpy.ndarray, shape=(3, 3), dtype=numpy.float32
        the box_vectors of the configuration. The coordinates are wrapped in a
        simtk.unit.Unit.
    """
    if snapshot.statics is not None:
        return unmask_quantity(snapshot.statics.box_vectors)

    return None


@box_vectors.setter
def box_vectors(self, value):
    if value is not None:
        sc = StaticContainer(box_vectors=value,
                             coordinates=self.coordinates,
                             engine=self.engine)
    else:
        sc = None

    self.statics = sc


@property
def md(snapshot):
    """
    Returns
    -------
    md : mdtraj.Trajectory
        the actual trajectory object. Can be used with all functions from mdtraj

    Notes
    -----
    Rather slow since the topology has to be made each time. Try to avoid
    it. This will only work if the engine has an mdtraj_topology property.
    """
    if snapshot.statics is not None:
        return paths.Trajectory([snapshot]).to_mdtraj()


@property
def xyz(snapshot):
    """
    Returns
    -------
    xyz : numpy.ndarray, shape=(atoms, 3), dtype=numpy.float32
        atomic coordinates without dimensions. Be careful.

    """
    import simtk.unit as u

    coord = snapshot.coordinates
    if type(coord) is u.Quantity:
        return coord._value
    else:
        return coord
