import numpy as np

@property
def n_degrees_of_freedom(snapshot):
    """
    Returns
    -------
    n_degrees_of_freedom: int
        number of degrees of freedom in this system (after accounting for
        constraints)
    """
    return snapshot.engine.n_degrees_of_freedom()

@property
def instantaneous_temperature(snapshot):
    """
    Returns
    -------
    instantaneous_temperature : float
        instantaneous temperature from the kinetic energy of this snapshot
    """
    masses = snapshot.masses
    velocities = snapshot.velocities
    double_ke = sum([masses[i] * velocities[i].dot(velocities[i])
                     for i in range(len(masses))])

    dofs = snapshot.n_degrees_of_freedom

    temperature = double_ke / dofs
    return temperature


