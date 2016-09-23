import numpy as np

@property
def n_degrees_of_freedom(snapshot):
    topol = snapshot.engine.topology
    return topol.n_atoms * topol.n_spatial

@property
def instantaneous_temperature(snapshot):
    masses = snapshot.masses
    velocities = snapshot.velocities
    double_ke = sum([masses[i] * velocities[i].dot(velocities[i])
                     for i in range(len(masses))])

    dofs = snapshot.n_degrees_of_freedom

    temperature = double_ke / dofs
    return temperature


