import os
import numpy as np
from integrators import VVVRIntegrator

import simtk.unit as u
import simtk.openmm as openmm
from simtk.openmm.app import ForceField, PME, HBonds, PDBFile

from opentis.storage import Storage
from opentis.tools import snapshot_from_pdb, to_openmm_topology
from opentis.wrapper import storable
from opentis.dynamics_engine import DynamicsEngine
from opentis.snapshot import Snapshot, Configuration, Momentum

@storable
class OpenMMEngine(DynamicsEngine):
    """We only need a few things from the simulation. This object duck-types
    an OpenMM simulation object so that it quacks the methods we need to
    use."""

    units = {
        'length' : u.nanometers,
        'velocity' : u.nanometers / u.picoseconds,
        'energy' : u.kilojoules_per_mole
    }

    @staticmethod
    def auto(filename, template, options, mode='auto'):
        if mode == 'auto':
            if os.path.isfile(filename):
                mode = 'restore'
            else:
                mode = 'create'

        if mode == 'create':
            return OpenMMEngine.create_with_storage(filename, template, options)
        elif mode == 'restore':
            return OpenMMEngine.restore_from_storage(filename)
        else:
            raise ValueError('Unknown mode: ' + mode)
            return None


    @staticmethod
    def create_with_storage(filename, template, options):
        # for openmm we will create a suitable for configuration with
        # attached box_vectors and topolgy

        if type(template) is str:
            template = snapshot_from_pdb(template)

        # once we have a template configuration (coordinates to not really matter)
        # we can create a storage. We might move this logic out of the dynamics engine
        # and keep storage and engine generation completely separate!

        storage = Storage(
            filename=filename,
            template=template,
            mode='w'
        )

        # save simulator options, should be replaced by just saving the simulator object
        storage.init_str('simulation_options')
        storage.write_as_json('simulation_options', options)

        engine = OpenMMEngine(template, options)
        engine.storage = storage

        return engine

    @staticmethod
    def restore_from_storage(filename):
        # open storage, which also gets the topology!
        storage = Storage(
            filename=filename,
            mode='a'
        )

        options = storage.restore_object('simulation_options')
        engine = OpenMMEngine(storage.template, options)

        engine.storage = storage

        return engine


    def __init__(self, template, options=None):

        if options is not None:
            self.options = options
        else:
            self.options = {}
        super(OpenMMEngine, self).__init__(
            options=self.options
        )

        self.template = template
        self.topology = template.topology

        # set up the max_length_stopper (if n_frames_max is given)
        if 'nsteps_per_frame' in self.options:
            self.nsteps_per_frame = self.options['nsteps_per_frame']

        if 'solute_indices' in self.options:
            self.solute_indices = self.options['solute_indices']

        if 'n_frames_max' in self.options:
            self.n_frames_max = self.options['n_frames_max']

        self.n_atoms = self.topology.n_atoms

        # set up the OpenMM simulation
        forcefield = ForceField( self.options["forcefield_solute"],
                                 self.options["forcefield_solvent"] )

        openmm_topology = to_openmm_topology(self.template)

        system = forcefield.createSystem( openmm_topology,
                                          nonbondedMethod=PME,
                                          nonbondedCutoff=1.0 * u.nanometers,
                                          constraints=HBonds )

#        self.system_serial = openmm.XmlSerializer.serialize(system)

        integrator = VVVRIntegrator( self.options["temperature"],
                                     self.options["collision_rate"],
                                     self.options["timestep"] )

#        self.integrator_serial = openmm.XmlSerializer.serialize(system)
#        self.integrator_class = type(integrator).__name__

        simulation = openmm.app.Simulation(openmm_topology, system,
                                           integrator)

        # claim the OpenMM simulation as our own
        self.simulation = simulation


    def equilibrate(self, nsteps):
        # TODO: rename... this is position restrained equil, right?
        #self.simulation.context.setPositions(self.pdb.positions) #TODO move
        system = self.simulation.system
        n_solute = len(self.solute_indices)

        solute_masses = u.Quantity(np.zeros(n_solute, np.double), u.dalton)
        for i in self.solute_indices:
            solute_masses[i] = system.getParticleMass(i)
            system.setParticleMass(i,0.0)

        self.simulation.step(nsteps)

        for i in self.solute_indices:
            system.setParticleMass(i, solute_masses[i].value_in_unit(u.dalton))



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