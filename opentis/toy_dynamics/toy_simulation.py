
class ToySimulation(object):
    '''The trick is that we have various "simulation" classes (either
    generated directly as here, or subclassed for more complication
    simulation objects as in OpenMM), but they all quack the same when it
    comes to things the Simulator calls on them for'''

    def get_velocities(self):
        pass

    def get_configurations(self):
        pass

    def get_snapshot(self):
        pass

    def get_kinetic_energy(self):
        pass

    def get_potential_energy(self):
        pass

    def get_periodic_box(self):
        pass

    def init_trajectory_with_snapshot(self, snapshot):
        pass

    def step(self, nsteps):
        pass

# TODO: mostly need to abstract out things in Simulator.py and snapshot.py
# into an OpenMMSimulation object that also looks like this one
