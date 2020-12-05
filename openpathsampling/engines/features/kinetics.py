from .shared import KineticContainerStore, KineticContainer, unmask_quantity
from openpathsampling.netcdfplus import WeakLRUCache

variables = ['kinetics', 'is_reversed']
lazy = ['kinetics']
flip = ['is_reversed']

dimensions = ['n_atoms', 'n_spatial']

_vel_unit = "simtk(unit.nanometer/unit.picosecond)"
_vel_str = "ndarray.float32({n_atoms},{n_spatial})"
schema_entries = [
    ('kinetics', [
        ('velocities', _vel_unit + "*" + _vel_str),
        ('engine', 'uuid'),
    ]),
    ('is_reversed', 'bool'),
]


def netcdfplus_init(store):
    kinetic_store = KineticContainerStore()
    kinetic_store.set_caching(WeakLRUCache(10000))

    name = store.prefix + 'kinetics'

    # tell the KineticContainerStore to base its dimensions on names
    # prefixed with the store name to make sure each snapshot store has
    # its own kinetics store
    kinetic_store.set_dimension_prefix_store(store)

    store.storage.create_store(name, kinetic_store, False)
    store.create_variable(
        'kinetics', 'lazyobj.' + name,
        description="the snapshot index (0..n_momentum-1) 'frame' of "
                    "snapshot '{idx}'.")

    store.create_variable(
        'is_reversed', 'bool',
        description="the indicator if momenta should be flipped.")


@property
def velocities(self):
    """
    The velocities in the configuration. If the snapshot is reversed a
    copy of the original (unreversed) velocities is made which is then
    returned
    """
    if self.kinetics is not None:
        vel = unmask_quantity(self.kinetics.velocities)
        if self.is_reversed:
            return -1.0 * vel
        else:
            return vel

    return None


@velocities.setter
def velocities(self, value):
    if value is not None:
        kc = KineticContainer(velocities=value, engine=self.engine)
    else:
        kc = None

    self.is_reversed = False
    self.kinetics = kc
