from shared import KineticContainerStore

variables = ['kinetics', 'is_reversed']
lazy = ['kinetics']
flip = ['is_reversed']


def netcdfplus_init(store):
    store.storage.create_store('kinetics', KineticContainerStore())

    store.create_variable('kinetics', 'lazyobj.kinetics',
                        description="the snapshot index (0..n_momentum-1) 'frame' of snapshot '{idx}'.",
                        )

    store.create_variable('is_reversed', 'bool',
                        description="the indicator if momenta should be flipped.",
                        )

@property
def velocities(self):
    """
    The velocities in the configuration. If the snapshot is reversed a
    copy of the original (unreversed) velocities is made which is then
    returned
    """
    if self.kinetics is not None:
        if self.is_reversed:
            return -1.0 * self.kinetics.velocities
        else:
            return self.kinetics.velocities

    return None
