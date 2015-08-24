'''
Created on 06.07.2014

@author: JDC Chodera
@author: JH Prinz
'''

import logging

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')

import openpathsampling as paths
import simtk.unit as u
from netcdfplus import NetCDFPlus

#=============================================================================================
# OPS SPECIFIC STORAGE
#=============================================================================================

class Storage(NetCDFPlus):

    def get_unit(self, dimension):
        """
        Return a simtk.Unit instance from the unit_system the is of the specified dimension, e.g. length, time
        """
        return u.Unit({self.unit_system.base_units[u.BaseDimension(dimension)] : 1.0})

    @property
    def template(self):
        """
        Return the template snapshot from the storage

        Returns
        -------
        Snapshot
            the initial snapshot
        """
        if self._template is None:
            self._template = self.snapshots.load(int(self.variables['template_idx'][0]))

        return self._template

    def _setup_class(self):
        super(Storage, self)._setup_class()
                # use MD units

        self.dimension_units = {
            'length': u.nanometers,
            'velocity': u.nanometers / u.picoseconds,
            'energy': u.kilojoules_per_mole
        }


    def clone(self, filename, subset):
        """
        Creates a copy of the netCDF file and allows to reduce the used atoms.

        Notes
        -----
        This is mostly used to remove water but keep the data intact.
        """

        storage2 = Storage(filename=filename, template=self.template.subset(subset), mode='w')

        # Copy all configurations and momenta to new file in reduced form
        for obj in self.configurations:
            storage2.configurations.save(obj.copy(subset), idx=obj.idx[self.configurations])
        for obj in self.momenta.iterator():
            storage2.momenta.save(obj.copy(subset), idx=obj.idx[self.momenta])

        # All other should be copied one to one. We do this explicitely although we could just copy all
        # and exclude configurations and momenta, but this seems cleaner

        for storage_name in [
                'pathmovers', 'topologies', 'networks', 'details', 'trajectories',
                'shootingpoints', 'shootingpointselectors', 'engines', 'volumes',
                'samplesets', 'ensembles', 'transitions', 'steps', 'pathmovechanges',
                'samples', 'snapshots', 'pathsimulators', 'cvs'
            ]:
            self.clone_storage(storage_name, storage2)

        storage2.close()

    # TODO: Need to copy cvs without caches!
    def clone_empty(self, filename):
        """
        Creates a copy of the netCDF file and replicates only the static parts which I consider
            ensembles, volumes, engines, pathmovers, shootingpointselectors. We do not need to
            reconstruct collectivevariables since these need to be created again completely and then
            the necessary arrays in the file will be created automatically anyway.

        Notes
        -----
        This is mostly used to restart with a fresh file. Same setup, no results.
        """
        storage2 = Storage(filename=filename, template=self.template, mode='w')

        for storage_name in [
                'pathmovers', 'topologies', 'networks',
                'shootingpointselectors', 'engines', 'volumes',
                'ensembles', 'transitions', 'pathsimulators'
            ]:
            self.clone_storage(storage_name, storage2)

        storage2.close()

    @property
    def n_atoms(self):
        return self.topology.n_atoms

    @property
    def n_spatial(self):
        return self.topology.n_spatial

    def __init__(self, filename, mode=None,
                 template=None, units=None):
        '''
        Create a storage for complex objects in a netCDF file

        Parameters
        ----------
        filename : string
            filename of the netcdf file to be used or created
        mode : string, default: None
            the mode of file creation, one of 'w' (write), 'a' (append) or
            None, which will append any existing files.
        template : openpathsampling.Snapshot
            a Snapshot instance that contains a reference to a Topology, the
            number of atoms and used units
        units : dict of {str : simtk.unit.Unit } or None
            representing a dict of string representing a dimension
            ('length', 'velocity', 'energy') pointing to
            the simtk.unit.Unit to be used. If not None overrides the
            standard units used
        '''

        self._template = template
        super(Storage, self).__init__(filename, mode, units=units)

    def _register_storages(self):
        """
        Register all Stores used in the OpenPathSampling Storage

        """

        # objects with special storages

        # self.objectname = ... could also be done in the initialization
        # automatically. But the IDE would not be able to autocomplete
        # so we leave it this way :)

        self.add('trajectories', paths.storage.TrajectoryStore())
        self.add('snapshots', paths.storage.SnapshotStore())
        self.add('configurations', paths.storage.ConfigurationStore())
        self.add('momenta', paths.storage.MomentumStore())
        self.add('samples', paths.storage.SampleStore())
        self.add('samplesets', paths.storage.SampleSetStore())
        self.add('pathmovechanges', paths.storage.PathMoveChangeStore())
        self.add('steps', paths.storage.MCStepStore())

        self.add('cvs', paths.storage.ObjectDictStore(paths.CollectiveVariable, paths.Snapshot))
        self.collectivevariables = self.cvs

        # normal objects

        self.add('details', paths.storage.ObjectStore(paths.Details, has_uid=False, has_name=False))
        self.add('topologies', paths.storage.ObjectStore(paths.Topology, has_uid=True, has_name=True))
        self.add('pathmovers', paths.storage.ObjectStore(paths.PathMover, has_uid=True, has_name=True))
        self.add('shootingpoints', paths.storage.ObjectStore(paths.ShootingPoint, has_uid=False))
        self.add('shootingpointselectors', paths.storage.ObjectStore(paths.ShootingPointSelector, has_uid=False, has_name=True))
        self.add('engines', paths.storage.ObjectStore(paths.DynamicsEngine, has_uid=True, has_name=True))
        self.add('pathsimulators', paths.storage.ObjectStore(paths.PathSimulator, has_uid=True, has_name=True))
        self.add('transitions', paths.storage.ObjectStore(paths.Transition, has_uid=True, has_name=True))
        self.add('networks', paths.storage.ObjectStore(paths.TransitionNetwork, has_uid=True, has_name=True))

        # nestable objects

        self.add('volumes', paths.storage.ObjectStore(paths.Volume, has_uid=True, nestable=True, has_name=True))
        self.add('ensembles', paths.storage.ObjectStore(paths.Ensemble, has_uid=True, nestable=True, has_name=True))

    def _initialize(self):
        # Set global attributes.
        setattr(self, 'title', 'OpenPathSampling Storage')
        setattr(self, 'ConventionVersion', '0.2')

        template = self._template

        if template.topology is not None:
            self.topology = template.topology
        else:
            raise RuntimeError("A Storage needs a template snapshot with a topology")

        if 'atom' not in self.dimensions:
            self.createDimension('atom', self.topology.n_atoms)

        # spatial dimensions
        if 'spatial' not in self.dimensions:
            self.createDimension('spatial', self.n_spatial)

        # update the units for dimensions from the template
        self.dimension_units.update(paths.tools.units_from_snapshot(template))
        self._init_storages()

        logger.info("Saving topology")

        topology_id = self.topologies.save(self.topology)

        logger.info("Create initial template snapshot")

        # Save the initial configuration
        self.snapshots.save(template)

        self.createVariable('template_idx', 'i4', 'scalar')
        self.variables['template_idx'][:] = template.idx[self.snapshots]

    def _restore(self):
        self.topology = self.topologies[0]


class AnalysisStorage(Storage):
    def __init__(self, filename):
        super(AnalysisStorage, self).__init__(
            filename=filename,
            mode='r'
        )

        self.set_caching_mode('analysis')

        # Let's go caching

        self.samples.cache_all()
        self.samplesets.cache_all()
        self.cvs.cache_all()
        map(lambda x : x.cache_all(self.cvs), self.cvs)
        self.volumes.cache_all()
        self.ensembles.cache_all()
        self.pathmovers.cache_all()
        self.pathmovechanges.cache_all()
        self.steps.cache_all()
