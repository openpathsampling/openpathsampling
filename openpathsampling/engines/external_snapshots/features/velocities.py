variables = ['velocity_direction']
minus = ['velocity_direction']
default_none = ['_velocities']

def netcdfplus_init(store):
    store.create_variable(
        'velocity_direction',
        'int',
        description=("whether the stored velocities are flipped. "
                     + "Note: all trajectories (forward or backward in "
                     + "in trajectory time) are generated as forward "
                     + "trajectories.")
    )

@property
def velocities(snapshot):
    """
    """
    if snapshot._velocities is None:
        snapshot.load_details()

    return snapshot.velocity_direction * snapshot._velocities
