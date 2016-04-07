import logging

import simtk.openmm
import simtk.unit as u
from simtk.openmm.app import Simulation

from openpathsampling.engines import DynamicsEngine
from snapshot import Snapshot

logger = logging.getLogger(__name__)


class OpenMMEngine(DynamicsEngine):
    """OpenMM dynamics engine based on using an 'simtk.openmm` system and integrator object.

    The engine will create a :class:`simtk.openmm.app.Simulation` instance and uses this to generate new frames.

    """

    units = {
        'length': u.nanometers,
        'velocity': u.nanometers / u.picoseconds,
        'energy': u.joule / u.mole
    }

    _default_options = {
        'nsteps_per_frame': 10,
        'n_frames_max': 5000,
        'platform': 'fastest'
    }

    base_snapshot_type = Snapshot

    #TODO: Planned to move topology to be part of engine and not snapshot
    #TODO: Deal with cases where we load a GPU based engine, but the platform is not available
    def __init__(self, template, system, integrator, options=None, properties=None):
        """
        Parameters
        ----------
        template : openpathsampling.Snapshot
            a template snapshots which provides the topology object to be used to create the openmm engine
        system : simtk.openmm.app.System
            the openmm system object
        integrator : simtk.openmm.Integrator
            the openmm integrator object
        options : dict
            a dictionary that provides additional settings for the OPS engine. Allowed are
                'n_steps_per_frame' : int, default: 10, the number of integration steps per returned snapshot
                'n_frames_max' : int or None, default: 5000, the maximal number of frames allowed for a returned
                trajectory object
                `platform` : str, default: `fastest`, the openmm specification for the platform to be used, also 'fastest' is allowed
                which will pick the currently fastest one available

        Notes
        -----
        the `n_frames_max` does not limit Trajectory objects in length. It only limits the maximal lenght of returned
        trajectory objects when this engine is used.
        picking `fasted` as platform will not save `fastest` as the platform but rather replace the platform with the
        currently fastest one (usually `OpenCL` or `CUDA` for GPU and `CPU` otherwise). If you load this engine it will
        assume the same engine and not the currently fastest one, so you might have to create a replacement that uses
        another engine.
        """

        self.system = system
        self.integrator = integrator

        super(OpenMMEngine, self).__init__(
            options=options,
            template=template
        )

        if self.options['platform'] == 'fastest':

            speed = 0.0
            platform = None

            # determine the fastest platform
            for platform_idx in range(simtk.openmm.Platform.getNumPlatforms()):
                pf = simtk.openmm.Platform.getPlatform(platform_idx)
                if pf.getSpeed() > speed:
                    speed = pf.getSpeed()
                    platform = pf.getName()

            if platform is not None:
                self.options['platform'] = platform

        if properties is None:
            properties = dict()

        self.properties = properties

        # set no cached snapshot, means it will be constructed from the openmm context
        self._current_snapshot = None
        self._current_momentum = None
        self._current_configuration = None
        self._current_box_vectors = None

        self._simulation = None

    def from_new_options(self, integrator=None, options=None):
        """
        Create a new engine with the same system, but different options and/or integrator

        Notes
        -----
        This can be used to quickly set up simulations at various temperatures or change the
        step sizes, etc...

        """
        if integrator is None:
            integrator = self.integrator

        new_options = dict()
        new_options.update(self.options)

        if options is not None:
            new_options.update(options)

        new_engine = OpenMMEngine(self.template, self.system, integrator, new_options)

        if integrator is self.integrator and new_engine.options['platform'] == self.options['platform']:
            # apparently we use a simulation object which is the same as the new one
            # since we do not change the platform or change the integrator
            # it means if it exists we copy the simulation object

            new_engine._simulation = self._simulation

        return new_engine

    @property
    def simulation(self):
        if self._simulation is None:
            self.initialize()

        return self._simulation

    def initialize(self):
        """
        Create the final OpenMMEngine

        Notes
        -----
        This step is OpenMM specific and will actually create the openmm.Simulation object used
        to run the simulations. The object will be created automatically the first time the
        engine is used. This way we will not create unnecessay Engines in memory during analysis.

        """

        if self._simulation is None:
            self._simulation = simtk.openmm.app.Simulation(
                topology=self.template.topology.md.to_openmm(),
                system=self.system,
                integrator=self.integrator,
                platform=simtk.openmm.Platform.getPlatformByName(self.platform)
            )

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
            'options' : self.options,
            'properties' : self.properties
        }

    @classmethod
    def from_dict(cls, dct):
        system_xml = dct['system_xml']
        integrator_xml = dct['integrator_xml']
        template = dct['template']
        options = dct['options']
        properties = dct['properties']

        return OpenMMEngine(
            template=template,
            system=simtk.openmm.XmlSerializer.deserialize(system_xml),
            integrator=simtk.openmm.XmlSerializer.deserialize(integrator_xml),
            options=options,
            properties=properties
        )

    @property
    def snapshot_timestep(self):
        return self.nsteps_per_frame * self.simulation.integrator.getStepSize()

    def _build_current_snapshot(self):
        # TODO: Add caching for this and mark if changed

        state = self.simulation.context.getState(getPositions=True,
                                                 getVelocities=True,
                                                 getEnergy=True)

        snapshot = Snapshot.construct(
            coordinates=state.getPositions(asNumpy=True),
            box_vectors=state.getPeriodicBoxVectors(asNumpy=True),
            velocities=state.getVelocities(asNumpy=True),
            topology=self.topology
        )

        return snapshot

    @property
    def current_snapshot(self):
        if self._current_snapshot is None:
            self._current_snapshot = self._build_current_snapshot()

        return self._current_snapshot

    def _changed(self):
        self._current_snapshot = None

    @current_snapshot.setter
    def current_snapshot(self, snapshot):

        self.check_snapshot_type(snapshot)

        if snapshot is not self._current_snapshot:
            if snapshot.coordinates is not None:
                self.simulation.context.setPositions(snapshot.coordinates)

            if snapshot.velocities is not None:
                self.simulation.context.setVelocities(snapshot.velocities)

            # After the updates cache the new snapshot
            self._current_snapshot = snapshot

    def generate_next_frame(self):
        self.simulation.step(self.nsteps_per_frame)
        self._current_snapshot = None
        return self.current_snapshot

    def minimize(self):
        self.simulation.minimizeEnergy()
        # make sure that we get the minimized structure on request
        self._current_snapshot = None

