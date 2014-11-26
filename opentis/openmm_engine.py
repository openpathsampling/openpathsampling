import numpy as np
from snapshot import Snapshot, Configuration, Momentum
from dynamics_engine import DynamicsEngine
from integrators import VVVRIntegrator

from simtk.unit import femtoseconds, picoseconds, nanometers, kelvin, dalton
from simtk.unit import Quantity

import simtk.unit as units

import simtk.openmm.app

import simtk.openmm as openmm
from simtk.openmm.app import ForceField, PME, HBonds, PDBFile
import mdtraj as md

from opentis.storage import Storage

from opentis.tools import snapshot_from_pdb, to_openmm_topology

import os

class OpenMMEngine(DynamicsEngine):
    """We only need a few things from the simulation. This object duck-types
    an OpenMM simulation object so that it quacks the methods we need to
    use."""

    def __init__(self, filename, topology_file, options, mode='auto'):

        if mode == 'auto':
            if os.path.isfile(filename):
                mode = 'restore'
            else:
                mode = 'create'

        self.fn_storage = filename

        if mode == 'create':

            self.options = options

            # for openmm we will create a suitable for configuration with
            # attached box_vectors and topolgy

            self.initial_configuration = None

            if isinstance(topology_file, md.Topology):
                self.topology = topology_file
                # This also needs box_dimensions
                # Right now this is not supported

            elif isinstance(topology_file, simtk.openmm.app.Topology):
                # In case of an openmm toppology we use the box size in the topologgy and create empty coordinates
                n_atoms = self.topology_file.n_atoms

                self.initial_configuration = Configuration(
                    coordinates=units.Quantity(np.zeros((n_atoms, 3)), units.nanometer),
                    box_vectors=units.Quantity(self.topology.setUnitCellDimensions(), units.nanometer),
                    potential_energy=units.Quantity(0.0, units.kilojoules_per_mole),
                    topology=md.Topology.from_openmm(topology_file)
                )

            elif type(topology_file) is str:
                self.initial_configuration = snapshot_from_pdb(topology_file)

            # once we have a template configuration (coordinates to not really matter)
            # we can create a storage. We might move this logic out of the dynamics engine
            # and keep sotrage and engine generation completely separate!

            storage = Storage(
                filename=self.fn_storage,
                template=self.initial_configuration,
                mode='w'
            )

            # save simulator options, should be replaced by just saving the simulator object
            storage.init_str('simulation_options')
            storage.write_as_json('simulation_options', self.options)

            self.storage = storage

        elif mode == 'restore':
            # open storage, which also gets the topology!
            self.storage = Storage(
                filename=self.fn_storage,
                mode='a'
            )

            self.options = self.storage.restore_object('simulation_options')

        # finally restored storage, topology and the options to lets build the actual engine

        self.topology = self.storage.topology
        self.initial_configuration = self.storage.initial_configuration

        super(OpenMMEngine, self).__init__(
            options=options
        )

        # set up the max_length_stopper (if n_frames_max is given)
        if 'nsteps_per_frame' in self.options:
            self.nsteps_per_frame = self.options['nsteps_per_frame']

        if 'solute_indices' in self.options:
            self.solute_indices = self.options['solute_indices']

        if 'n_frames_max' in self.options:
            self.n_frames_max = self.options['n_frames_max']

        self.n_atoms = self.topology.n_atoms

        # set up the OpenMM simulation
        forcefield = ForceField( options["forcefield_solute"],
                                 options["forcefield_solvent"] )

        openmm_topology = to_openmm_topology(self.initial_configuration)

        system = forcefield.createSystem( openmm_topology,
                                          nonbondedMethod=PME,
                                          nonbondedCutoff=1.0 * nanometers,
                                          constraints=HBonds )

#        self.system_serial = openmm.XmlSerializer.serialize(system)

        integrator = VVVRIntegrator( options["temperature"],
                                     options["collision_rate"],
                                     options["timestep"] )

#        self.integrator_serial = openmm.XmlSerializer.serialize(system)
#        self.integrator_class = type(integrator).__name__

        simulation = openmm.app.Simulation(openmm_topology, system,
                                           integrator)

        # claim the OpenMM simulation as our own
        self.simulation = simulation

        return


    def equilibrate(self, nsteps):
        # TODO: rename... this is position restrained equil, right?
        #self.simulation.context.setPositions(self.pdb.positions) #TODO move
        system = self.simulation.system
        n_solute = len(self.solute_indices)

        solute_masses = Quantity(np.zeros(n_solute, np.double), dalton)
        for i in self.solute_indices:
            solute_masses[i] = system.getParticleMass(i)
            system.setParticleMass(i,0.0)

        self.simulation.step(nsteps)

        for i in self.solute_indices:
            system.setParticleMass(i, solute_masses[i].value_in_unit(dalton))



    # this property is specific to direct control simulations: other
    # simulations might not use this
    # TODO: Maybe remove this and put it into the creation logic
    @property
    def nsteps_per_frame(self):
        return self._nsteps_per_frame

    @nsteps_per_frame.setter
    def nsteps_per_frame(self, value):
        self._nsteps_per_frame = value

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
                        box_vectors = state.getPeriodicBoxVectors(asNumpy=True),
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
        self.configuration = snapshot.configuration
        self.momentum = snapshot.momentum

    def generate_next_frame(self):
        self.simulation.step(self.nsteps_per_frame)
        return self.current_snapshot


    # (possibly temporary) shortcuts for momentum and configuration
    @property
    def momentum(self):
        state = self.simulation.context.getState(getVelocities=True,
                                                 getEnergy=True)
        return Momentum(velocities = state.getVelocities(asNumpy=True),
                        kinetic_energy = state.getKineticEnergy()
                       )

    @momentum.setter
    def momentum(self, momentum):
        self.simulation.context.setVelocities(momentum.velocities)

    @property
    def configuration(self):
        state = self.simulation.context.getState(getPositions=True,
                                                 getVelocities=True,
                                                 getEnergy=True)
        return Configuration(coordinates = state.getPositions(asNumpy=True),
                             box_vectors = state.getPeriodicBoxVectors(asNumpy=True),
                             potential_energy = state.getPotentialEnergy(),
                             topology = self.topology
                            )

    @configuration.setter
    def configuration(self, config):
        self.simulation.context.setPositions(config.coordinates)
        # TODO: Check if this is the right way to make sure the box is right!
        # self.simulation.context.getPeriodicBoxVectors(config.box_vectors)