
from simtk.openmm.app import Simulation
from snapshot import Snapshot

class OpenMMSimulation(Simulation):
    """We only need a few things from the simulation. This object duck-types
    an OpenMM simulation object so that it quacks the methods we need to
    use."""

    # this property is specific to certain simulations: other simulations
    # might not use this
    @property
    def nsteps_per_iteration(self):
        return self._nsteps_per_iteration

    @nsteps_per_iteration.setter
    def nsteps_per_iteration(self, value):
        self._nsteps_per_iteration = value


    def load_momentum(self, momentum):
        """Loads the current context velocities/kinetic energy into the
        `Momentum` object given as a parameter

        Parameters
        ----------
        momentum : Momentum
            output object into which the current simulation state is loaded
        """
        state = self.context.getState(getVelocities=True, getEnergy=True)
        momentum._velocities = state.getVelocities(asNumpy=True)
        momentum._kinetic_energy = state.getKineticEnergy()

    def load_configuration(self, configuration):
        """Loads the current context positions/potential energy/box vectors
        into the `Configuration` object given as a parameter

        Parameters
        ----------
        configuration : Configuration
            output object into which the current simulation state is loaded
        """
        # this is basically the stuff that used to be in
        # Configuration.__init__ in the case that a context was given
        state = self.context.getState(getPositions=True, getEnergy=True)
        configuration._coordinates = state.getPositions(asNumpy=True)
        configuration._box_vectors = state.getPeriodicBoxVectors()
        configuration._potential_energy = state.getPotentialEnergy()

    def load_snapshot(self, snapshot):
        """Loads the current context into the `Snapshot` object given as a
        parameter. Include positions, velocities, kinetic and potential
        energy, and box vectors.

        Parameters
        ----------
        snapshot : Snapshot
            output object into which the current simulation state is loaded
        """
        state = self.context.getState(getPositions=True, getVelocities=True,
                                      getEnergy=True)
        snapshot.configuration._coordinates = state.getPositions(asNumpy=True)
        snapshot.configuration._box_vectors = state.getPeriodicBoxVectors()
        snapshot.configuration._potential_energy = state.getPotentialEnergy()
        snapshot.momentum._velocities = state.getVelocities(asNumpy=True)
        snapshot.momentum._kinetic_energy = state.getKineticEnergy()

    def init_simulation_with_snapshot(self, snapshot):
        self.context.setPositions(snapshot.coordinates)
        self.context.setVelocities(snapshot.velocities)


    def generate_next_frame(self):
        self.step(self.nsteps_per_iteration)
        return Snapshot(self)

    def stop(self):
        """Nothing special needs to be done to an OpenMMSimulation when you
        hit a stop condition."""
        pass
