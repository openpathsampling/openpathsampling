from shared import MomentumStore

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
