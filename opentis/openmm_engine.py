import numpy as np
import simtk.openmm
from snapshot import Snapshot, Configuration, Momentum
from trajectory import Trajectory
from dynamics_engine import DynamicsEngine
from ensemble import LengthEnsemble
from integrators import VVVRIntegrator
from storage import Storage

from simtk.unit import femtoseconds, picoseconds, nanometers, kelvin, dalton
from simtk.unit import Quantity

from simtk.openmm.app.pdbfile import PDBFile
import simtk.openmm as openmm
from simtk.openmm.app import ForceField, PME, HBonds

# TODO: figure out how much of this can be moved into the general case
class OpenMMEngine(DynamicsEngine):
    """We only need a few things from the simulation. This object duck-types
    an OpenMM simulation object so that it quacks the methods we need to
    use."""

    def __init__(self, filename, topology_file, opts, mode='auto'):
        super(OpenMMEngine, self).__init__()

        # tell everybody who their engine is
        Snapshot.engine = self
        Configuration.engine = self
        Trajectory.engine = self

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
            self.storage.engine = self
            self.storage._store_options(self)
            Trajectory.storage = self.storage
        if mode == 'restore':
            pass
        return


    def add_stored_parameters(self, param_dict):
        '''Adds parameters in param_dict to the attribute dictionary of the
        DynamicsEngine object, and saves the relevant keys as options to
        store.

        Parameters
        ----------
        param_dict : dict
            dictionary of attributes to be added (and stored); attribute
            names are keys, with appropriate values
        '''
        # TODO: I think this should go into the DynamicsEngine object
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

    # TODO: these two become defaults in DynamicsEngine() when 
    def start(self, snapshot=None):
        if snapshot is not None:
            self.current_snapshot = snapshot

    def stop(self, trajectory):
        """Nothing special needs to be done to an OpenMMSimulation when you
        hit a stop condition."""
        pass
