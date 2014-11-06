
class ToySimulation(object):
    '''The trick is that we have various "simulation" classes (either
    generated directly as here, or subclassed for more complication
    simulation objects as in OpenMM), but they all quack the same when it
    comes to things the Simulator calls on them for'''

    def __init__(self, pes, integ):
        self.pes = pes
        self.integ = integ

    @property
    def nsteps_per_iteration(self):
        return self._nsteps_per_iteration

    @nsteps_per_iteration.setter
    def nsteps_per_iteration(self, value):
        self._nsteps_per_iteration = value

    def load_momentum(self, momentum):
        momentum._velocities = self._velocities
        momentum._kinetic_energy = self.pes.kinetic_energy(self)

    def load_configuration(self, configuration):
        configuration._coordinates = self._coordinates
        configuration._potential_energy = self.pes.potential_energy(self)
        pass

    def load_snapshot(self, snapshot):
        pass

    def init_simulation_with_snapshot(self, snapshot):
        self.positions = snapshot.configuration.coordinates
        self.velocities = snap
        pass

    def generate_next_frame(self):
        pass

    def stop(self, trajectory):
        pass # pragma: no cover (no need to test this one)

