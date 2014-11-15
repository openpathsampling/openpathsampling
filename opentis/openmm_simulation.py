import numpy as np
import simtk.openmm
from snapshot import Snapshot, Configuration, Momentum
from trajectory import Trajectory
from Simulator import Simulator
from ensemble import LengthEnsemble
from integrators import VVVRIntegrator
from storage import Storage

from simtk.unit import femtoseconds, picoseconds, nanometers, kelvin, dalton
from simtk.unit import Quantity

from simtk.openmm.app.pdbfile import PDBFile
import simtk.openmm as openmm
from simtk.openmm.app import ForceField, PME, HBonds

# TODO: figure out how much of this can be moved into the general case
class OpenMMSimulation(Simulator):
    """We only need a few things from the simulation. This object duck-types
    an OpenMM simulation object so that it quacks the methods we need to
    use."""

    def __init__(self, filename, topology_file, opts, mode='auto'):
        super(OpenMMSimulation, self).__init__()

        # tell everybody who their simulator is
        Snapshot.simulator = self
        Configuration.simulator = self
        Trajectory.simulator = self

        # set up the opts
        self.opts = {}
        self.add_stored_parameters(opts)

        # storage
        self.fn_storage = filename 

        # topology
        self.pdb = PDBFile(topology_file)
        self.topology = self.pdb.topology

        if mode == 'create':
            # set up the OpenMM simulation
            forcefield = ForceField(self.forcefield_solute,
                                    self.forcefield_solvent)
            system = forcefield.createSystem(self.topology, 
                                             nonbondedMethod=PME,
                                             nonbondedCutoff=1*nanometers,
                                             constraints=HBonds)
            self.system_serial = openmm.XmlSerializer.serialize(system)

            integrator = VVVRIntegrator(self.temperature,
                                        self.collision_rate,
                                        self.timestep)
            self.integrator_serial = openmm.XmlSerializer.serialize(system)
            self.integrator_class = type(integrator).__name__

            simulation = openmm.app.Simulation(self.topology, system, 
                                               integrator)

            # claim the OpenMM simulation as our own
            self.simulation = simulation

            # set up the max_length_stopper (if n_frames_max is given)
            self.max_length_stopper = LengthEnsemble(slice(0,self.n_frames_max-1))

            # storage
            self.storage = Storage(
                topology_file=self.topology,
                filename=self.fn_storage,
                mode='w'
            )
            self.storage.simulator = self
            self.storage._store_options(self)
            Trajectory.storage = self.storage
        if mode == 'restore':
            pass
        return


    def add_stored_parameters(self, param_dict):
        '''Adds parameters in param_dict to the attribute dictionary of the
        simulator object, and saves the relevant keys as options to store.

        Parameters
        ----------
        param_dict : dict
            dictionary of attributes to be added (and stored); attribute
            names are keys, with appropriate values
        '''
        # TODO: I think this should go into the Simulator object
        for key in param_dict.keys():
            self.opts[key] = param_dict[key]
            setattr(self, key, param_dict[key])
        self.options_to_store = self.opts.keys()
        return
        

    def equilibrate(self, nsteps):
        # TODO: rename... this is position restrained equil, right?
        self.simulation.context.setPositions(self.pdb.positions)
        system = self.simulation.system
        n_solute = len(self.solute_indices)

        solute_masses = Quantity(np.zeros(n_solute, np.double), dalton)
        for i in self.solute_indices:
            solute_masses[i] = system.getParticleMass(i)
            system.setParticleMass(i,0.0)

        self.simulation.step(nsteps)

        for i in self.solute_indices:
            system.setParticleMass(i, solute_masses[i].value_in_unit(dalton))



    # this property is specific to certain simulations: other simulations
    # might not use this
    @property
    def nsteps_per_frame(self):
        return self._nsteps_per_frame

    @nsteps_per_frame.setter
    def nsteps_per_frame(self, value):
        self._nsteps_per_frame = value


    def load_momentum(self, momentum):
        """Loads the current context velocities/kinetic energy into the
        `Momentum` object given as a parameter

        Parameters
        ----------
        momentum : Momentum
            output object into which the current simulation state is loaded
        """
        state = self.simulation.context.getState(getVelocities=True,
                                                 getEnergy=True)
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
        state = self.simulation.context.getState(getPositions=True, 
                                                 getEnergy=True)
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
        state = self.simulation.context.getState(getPositions=True,
                                                 getVelocities=True,
                                                 getEnergy=True)
        snapshot.configuration._coordinates = state.getPositions(asNumpy=True)
        snapshot.configuration._box_vectors = state.getPeriodicBoxVectors()
        snapshot.configuration._potential_energy = state.getPotentialEnergy()
        snapshot.momentum._velocities = state.getVelocities(asNumpy=True)
        snapshot.momentum._kinetic_energy = state.getKineticEnergy()

    # TODO: there are two reasonable approaches to this: 
    # 1. require that part of the engine.next_frame() function be that the
    #    user saves a snapshot object called `self._current_snapshot`
    # 2. build the current snapshot on the fly every time the snapshot is
    #    needed
    # The trade-off is that (1) will be faster if we ask for the snapshot
    # frequently, but it is also much more likely to be a source of errors
    # for users who forget to implement that last step.

    def _build_current_snapshot(self):
        state = self.simulation.context.getState(getPositions=True,
                                                 getVelocities=True,
                                                 getEnergy=True)
        return Snapshot(coordinates = state.getPositions(asNumpy=True),
                        box_vectors = state.getPeriodicBoxVectors(),
                        potential_energy = state.getPotentialEnergy(),
                        velocities = state.getVelocities(asNumpy=True),
                        kinetic_energy = state.getKineticEnergy(),
                        topology = self.topology
                       )

    @property
    def current_snapshot(self):
        return self._build_current_snapshot()

    @current_snapshot.setter
    def current_snapshot(self, snapshot):
        self.simulation.context.setPositions(snapshot.coordinates)
        self.simulation.context.setVelocities(snapshot.velocities)

    def generate_next_frame(self):
        self.simulation.step(self.nsteps_per_frame)
        return self.current_snapshot

    # TODO: these two become defaults in Simulator() when 
    def start(self, snapshot=None):
        if snapshot is not None:
            self.current_snapshot = snapshot

    def stop(self, trajectory):
        """Nothing special needs to be done to an OpenMMSimulation when you
        hit a stop condition."""
        pass
