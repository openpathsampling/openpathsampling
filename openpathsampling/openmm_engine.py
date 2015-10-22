import os
import numpy as np
import simtk.unit as u
from simtk.openmm.app import ForceField, PME, HBonds, PDBFile, Simulation

import openpathsampling as paths
from openpathsampling.storage import Storage
from openpathsampling.integrators import VVVRIntegrator

class OpenMMRandomEngine(paths.DynamicsEngine):
    _default_options = {}

    def __init__(self, template=None):
        self.topology = template.topology
        self.options = {
        }

        super(OpenMMRandomEngine, self).__init__(
            options={},
            template=template
        )

        self.initialized = True

    def _build_current_snapshot(self):
        # TODO: Add caching for this and mark if changed

        tmp = self.template

        coordinates = u.Quantity(
            tmp.coordinates._value + np.random.normal(0.0, 0.02, tmp.coordinates.shape),
            tmp.coordinates.unit)
        velocities = u.Quantity(
            np.random.normal(0.0, 0.02, tmp.velocities.shape),
            tmp.velocities.unit)

        return paths.Snapshot(coordinates = coordinates,
                        box_vectors = tmp.box_vectors,
                        potential_energy = tmp.potential_energy,
                        velocities = velocities,
                        kinetic_energy = tmp.kinetic_energy,
                        topology = self.topology
                       )

    @property
    def current_snapshot(self):
        if self._current_snapshot is None:
            self._current_snapshot = self._build_current_snapshot()

        return self._current_snapshot

    @current_snapshot.setter
    def current_snapshot(self, snapshot):
        self._current_snapshot = snapshot

    def generate_next_frame(self):
        self._current_snapshot = None
        return self.current_snapshot


class OpenMMEngine(paths.DynamicsEngine):
    """OpenMM dynamics engine."""

    units = {
        'length': u.nanometers,
        'velocity': u.nanometers / u.picoseconds,
        'energy': u.joule / u.mole
    }

    _default_options = {
        'nsteps_per_frame': 10,
        'solute_indices': [0],
        'n_frames_max': 5000,
        "temperature": 300.0 * u.kelvin,
        'collision_rate': 1.0 / u.picoseconds,
        'timestep': 2.0 * u.femtoseconds,
        'platform': 'fastest',
        'forcefield_solute': 'amber96.xml',
        'forcefield_solvent': 'tip3p.xml'
    }

    @staticmethod
    def auto(filename, template, options, mode='auto', units=None):
        """
        Create or Restore a OpenMMEngine

        Parameters
        ----------
        filename : str
            the filename of the storage
        template : Snapshot
            the template Snapshot to be used. It contains the necessary
            units as well as the topology
        options : dict of { str : str }
            a dictionary that contains the parameters used in the
            construction of the OpenMMEngine
        mode : str ('restore', 'create' or 'auto')
            a string setting the mode of creation or restoration. The option
            'auto' (default) will only create a new storage if the file does
            not exist yet
        units : dict of {str : simtk.unit.Unit } or None (default)
            representing a dict of string representing a dimension
            ('length', 'velocity', 'energy') pointing the the
            simtk.unit.Unit to be used. This overrides the units used in the
            template

        Returns
        -------
        OpenMMEngine
            the created or restored OpenMMEngine instance

        """
        if mode == 'auto':
            if os.path.isfile(filename):
                mode = 'restore'
            else:
                mode = 'create'

        if mode == 'create':
            return OpenMMEngine._create_with_storage(filename, template, options, units)
        elif mode == 'restore':
            return OpenMMEngine._restore_from_storage(filename)
        else:
            raise ValueError('Unknown mode: ' + mode)
            return None


    @staticmethod
    def _create_with_storage(filename, template, options, units=None):
        # for openmm we will create a suitable for configuration with
        # attached box_vectors and topology

        if type(template) is str:
            template = paths.snapshot_from_pdb(template, units=units)

        # once we have a template configuration (coordinates to not really
        # matter) we can create a storage. We might move this logic out of
        # the dynamics engine and keep storage and engine generation
        # completely separate!

        storage = Storage(
            filename=filename,
            template=template,
            mode='w'
        )

        # save simulator options, should be replaced by just saving the simulator object

        options['template'] = template

        engine = OpenMMEngine(
            options=options
        )
        engine.storage = storage
        storage.engines.save(engine)

        return engine

    @staticmethod
    def _restore_from_storage(filename):
        # open storage, which also gets the topology!
        storage = Storage(
            filename=filename,
            mode='a'
        )

        engine = storage.engines.load(0)

        engine.storage = storage
        return engine


    def __init__(self, options, template=None):

        if 'template' in options:
            template = options['template']

        self.topology = template.topology
        self.options = {
        }

        super(OpenMMEngine, self).__init__(
            options=options,
            template=template
        )

        # set up the OpenMM simulation
        forcefield = ForceField( self.options["forcefield_solute"],
                                 self.options["forcefield_solvent"] )

        openmm_topology = paths.to_openmm_topology(self.template)

        system = forcefield.createSystem( openmm_topology,
                                          nonbondedMethod=PME,
                                          nonbondedCutoff=1.0 * u.nanometers,
                                          constraints=HBonds )

        integrator = VVVRIntegrator( self.options["temperature"],
                                     self.options["collision_rate"],
                                     self.options["timestep"] )

        simulation = Simulation(openmm_topology, system,
                                           integrator)

        # claim the OpenMM simulation as our own
        self.simulation = simulation

        # set no cached snapshot, menas it will be constructed from the openmm context
        self._current_snapshot = None
        self._current_momentum = None
        self._current_configuration = None
        self._current_box_vectors = None

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

        # empty cache
        self._current_snapshot = None

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

    @property
    def snapshot_timestep(self):
        return self.nsteps_per_frame * self.options['timestep']

    def _build_current_snapshot(self):
        # TODO: Add caching for this and mark if changed

        state = self.simulation.context.getState(getPositions=True,
                                                 getVelocities=True,
                                                 getEnergy=True)

        return paths.Snapshot(coordinates = state.getPositions(asNumpy=True),
                        box_vectors = state.getPeriodicBoxVectors(asNumpy=True),
                        potential_energy = state.getPotentialEnergy(),
                        velocities = state.getVelocities(asNumpy=True),
                        kinetic_energy = state.getKineticEnergy(),
                        topology = self.topology
                       )

    @property
    def current_snapshot(self):
        if self._current_snapshot is None:
            self._current_snapshot = self._build_current_snapshot()

        return self._current_snapshot

    def _changed(self):
        self._current_snapshot = None

    @current_snapshot.setter
    def current_snapshot(self, snapshot):
        if snapshot is not self._current_snapshot:
            if snapshot.configuration is not None:
                if self._current_snapshot is None or snapshot.configuration is not self._current_snapshot.configuration:
                    # new snapshot has a different configuration so update
                    self.simulation.context.setPositions(snapshot.coordinates)

                    # TODO: Check if this is the right way to make sure the box is right!
#                    if self._current_snapshot is None or snapshot.box_vectors != self._current_snapshot.box_vectors:
#                        self.simulation.context.getPeriodicBoxVectors(snapshot.box_vectors)

            if snapshot.momentum is not None:
                if self._current_snapshot is None or snapshot.momentum is not self._current_snapshot.momentum or snapshot.is_reversed != self._current_snapshot.is_reversed:
                    # new snapshot has a different momenta (different coordinates and reverse setting)
                    # so update. Note snapshot.velocities is different from snapshot.momenta.velocities!!!
                    # The first includes the reversal setting in the snapshot the second does not.
#                    print snapshot.momentum.velocities
                    self.simulation.context.setVelocities(snapshot.velocities)

            # After the updates cache the new snapshot
            self._current_snapshot = snapshot

    def generate_next_frame(self):
        self.simulation.step(self.nsteps_per_frame)
        self._current_snapshot = None
        return self.current_snapshot

    # (possibly temporary) shortcuts for momentum and configuration
    @property
    def momentum(self):
        return self.current_snapshot.momentum

    @property
    def configuration(self):
        return self.current_snapshot.configuration