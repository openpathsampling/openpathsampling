import numpy as np
from opentis.snapshot import Snapshot

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

    @property
    def mass(self):
        return self._mass
    @mass.setter
    def mass(self, value):
        self._mass = value
        self.minv = np.reciprocal(value)



    def load_momentum(self, momentum):
        momentum._velocities = self.velocities
        momentum._kinetic_energy = self.pes.kinetic_energy(self)

    def load_configuration(self, configuration):
        configuration._coordinates = self.positions
        configuration._potential_energy = self.pes.V(self)
        snapshot.configuration._box_vectors = None # toys without PBCs

    def load_snapshot(self, snapshot):
        snapshot.configuration._coordinates = self.positions
        snapshot.configuration._potential_energy = self.pes.V(self)
        snapshot.momentum._velocities = self.velocities
        snapshot.momentum._kinetic_energy = self.pes.kinetic_energy(self)
        snapshot.configuration._box_vectors = None

    def init_simulation_with_snapshot(self, snapshot):
        self.positions = snapshot.configuration.coordinates
        self.velocities = snapshot.momentum.velocities

    def generate_next_frame(self):
        self.integ.step(self, self.nsteps_per_iteration)
        return Snapshot(self)

    def stop(self, trajectory):
        pass # pragma: no cover (no need to test this one)

