import numpy as np
from snapshot import Snapshot, Configuration, Momentum
from trajectory import Trajectory
from dynamics_engine import DynamicsEngine
from ensemble import LengthEnsemble
from integrators import VVVRIntegrator
from storage import Storage

from simtk.unit import femtoseconds, picoseconds, nanometers, kelvin, dalton
from simtk.unit import Quantity

import simtk.openmm as openmm
from simtk.openmm.app import ForceField, PME, HBonds
import mdtraj as md

class OpenMMEngine(DynamicsEngine):
    """We only need a few things from the simulation. This object duck-types
    an OpenMM simulation object so that it quacks the methods we need to
    use."""

    def __init__(self, filename, topology_file, opts, mode='auto'):
        # if topology exists, it must be defined before running
        # super.__init__. This is ugly, but I don't see an easy way out.
        self.pdb = md.load(topology_file)
        topology = self.pdb.topology
        opts['topology'] = topology

        super(OpenMMEngine, self).__init__(filename=filename,
                                           options=opts,
                                           mode=mode)

        if mode == 'create':
            # set up the OpenMM simulation
            forcefield = ForceField( opts["forcefield_solute"],
                                     opts["forcefield_solvent"] )

            print 'CREATE'

            print type(self.pdb.topology)

            openmm_topology = topology.to_openmm()

            system = forcefield.createSystem( openmm_topology,
                                              nonbondedMethod=PME,
                                              nonbondedCutoff=1.0 * nanometers,
                                              constraints=HBonds )
            self.system_serial = openmm.XmlSerializer.serialize(system)

            integrator = VVVRIntegrator( opts["temperature"],
                                         opts["collision_rate"],
                                         opts["timestep"] )

            self.integrator_serial = openmm.XmlSerializer.serialize(system)
            self.integrator_class = type(integrator).__name__

            simulation = openmm.app.Simulation(openmm_topology, system,
                                               integrator)

            # claim the OpenMM simulation as our own
            self.simulation = simulation

        if mode == 'restore':
            pass
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
                             box_vectors = state.getPeriodicBoxVectors(),
                             potential_energy = state.getPotentialEnergy(),
                             topology = self.topology
                            )

    @configuration.setter
    def configuration(self, config):
        self.simulation.context.setPositions(config.coordinates)
