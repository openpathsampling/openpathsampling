import numpy as np
import simtk.unit as u
from simtk.openmm.app import ForceField, PME, HBonds, Simulation

import openpathsampling as paths
from openpathsampling.integrators import VVVRIntegrator

import simtk.openmm

class OpenMMEngine(paths.DynamicsEngine):
    """OpenMM dynamics engine based on a openmmtools.testsystem object.

    This is only to allow to use the examples from openmmtools.testsystems
    within this framework
    """

    units = {
        'length': u.nanometers,
        'velocity': u.nanometers / u.picoseconds,
        'energy': u.joule / u.mole
    }

    _default_options = {
        'nsteps_per_frame': 10,
        'n_frames_max': 5000,
        'platform' : 'CUDA'
    }

    def __init__(self, template, system, integrator, options=None):

        self.system = system
        self.integrator = integrator

        super(OpenMMEngine, self).__init__(
            options=options,
            template=template
        )

        # set no cached snapshot, means it will be constructed from the openmm context
        self._current_snapshot = None
        self._current_momentum = None
        self._current_configuration = None
        self._current_box_vectors = None

        self.simulation = None

    def create(self):
        """
        Create the final OpenMMEngine

        """
        self.simulation = simtk.openmm.app.Simulation(
            topology=self.template.topology.md.to_openmm(),
            system=self.system,
            integrator=self.integrator,
            platform=simtk.openmm.Platform.getPlatformByName(self.platform)
        )

        self.initialized = True


    @property
    def platform(self):
        return self.options['platform']

    def to_dict(self):
        system_xml = simtk.openmm.XmlSerializer.serialize(self.system)
        integrator_xml = simtk.openmm.XmlSerializer.serialize(self.integrator)

        return {
            'system_xml' : system_xml,
            'integrator_xml' : integrator_xml,
            'template' : self.template,
            'options' : self.options
        }

    @classmethod
    def from_dict(cls, dct):
        system_xml = dct['system_xml']
        integrator_xml = dct['integrator_xml']
        template = dct['template']
        options = dct['options']

        return cls(
            template=template,
            system=simtk.openmm.XmlSerializer.deserialize(system_xml),
            integrator=simtk.openmm.XmlSerializer.deserialize(integrator_xml),
            options=options
        )

    # this property is specific to direct control simulations: other
    # simulations might not use this
    # TODO: Maybe remove this and put it into the creation logic
    @property
    # this property is specific to direct control simulations: other
    # simulations might not use this
    # TODO: Maybe remove this and put it into the creation logic
    @property
    def nsteps_per_frame(self):
        return self._nsteps_per_frame

    @nsteps_per_frame.setter
    def nsteps_per_frame(self, value):
        self._nsteps_per_frame = value

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

            if snapshot.momentum is not None:
                if self._current_snapshot is None or snapshot.momentum is not self._current_snapshot.momentum or snapshot.is_reversed != self._current_snapshot.is_reversed:
                    self.simulation.context.setVelocities(snapshot.velocities)

            # After the updates cache the new snapshot
            self._current_snapshot = snapshot

    def generate_next_frame(self):
        self.simulation.step(self.nsteps_per_frame)
        self._current_snapshot = None
        return self.current_snapshot

    @property
    def momentum(self):
        return self.current_snapshot.momentum

    @property
    def configuration(self):
        return self.current_snapshot.configuration


class SimpleOpenMMEngine(OpenMMEngine):
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
        'platform': 'CUDA',
        'forcefield_solute': 'amber96.xml',
        'forcefield_solvent': 'tip3p.xml'
    }

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
                                          constraints=HBonds)

        integrator = VVVRIntegrator( self.options["temperature"],
                                     self.options["collision_rate"],
                                     self.options["timestep"])

        # claim the OpenMM simulation as our own
        self.system = system
        self.integrator = integrator

        # set no cached snapshot, menas it will be constructed from the openmm context
        self._current_snapshot = None
        self._current_momentum = None
        self._current_configuration = None
        self._current_box_vectors = None

        self.simulation = None

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