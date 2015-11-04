"""
Created on 06.07.2014

@author: JDC Chodera
@author: JH Prinz
"""

import logging

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')

import openpathsampling as paths
import simtk.unit as u
from netcdfplus import NetCDFPlus
from cache import WeakLRUCache, WeakValueCache


# =============================================================================================
# OPS SPECIFIC STORAGE
# =============================================================================================

class Storage(NetCDFPlus):
    """
    A netCDF4 wrapper to store trajectories based on snapshots of an OpenMM
    simulation. This allows effective storage of shooting trajectories
    """

    def get_unit(self, dimension):
        """
        Return a simtk.Unit instance from the unit_system the is of the specified dimension, e.g. length, time
        """
        return u.Unit({self.unit_system.base_units[u.BaseDimension(dimension)]: 1.0})

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
            storage2.configurations.save(obj.copy(subset=subset), idx=self.configurations.index[obj])
        for obj in self.momenta:
            storage2.momenta.save(obj.copy(subset=subset), idx=self.momenta.index[obj])

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
            ensembles, volumes, engines, path movers, shooting point selectors. We do not need to
            reconstruct collective variables since these need to be created again completely and then
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
        """
        Create a netdfplus storage for OPS Objects

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
        """

        self._template = template
        super(Storage, self).__init__(filename, mode, units=units)

    def _register_storages(self):
        """
        Register all Stores used in the OpenPathSampling Storage

        """

        # objects with special storages

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

        self.add('details', paths.storage.ObjectStore(paths.Details, has_name=False))
        self.add('topologies', paths.storage.ObjectStore(paths.Topology, has_name=True))
        self.add('pathmovers', paths.storage.ObjectStore(paths.PathMover, has_name=True))
        self.add('shootingpoints',
                 paths.storage.ObjectStore(paths.ShootingPoint, has_name=False))
        self.add('shootingpointselectors',
                 paths.storage.ObjectStore(paths.ShootingPointSelector, has_name=True))
        self.add('engines', paths.storage.ObjectStore(paths.DynamicsEngine, has_name=True))
        self.add('pathsimulators',
                 paths.storage.ObjectStore(paths.PathSimulator, has_name=True))
        self.add('transitions', paths.storage.ObjectStore(paths.Transition, has_name=True))
        self.add('networks',
                 paths.storage.ObjectStore(paths.TransitionNetwork, has_name=True))
        self.add('schemes',
                 paths.storage.ObjectStore(paths.MoveScheme, has_name=True))

        # nestable objects

        self.add('volumes',
                 paths.storage.ObjectStore(paths.Volume, nestable=True, has_name=True))
        self.add('ensembles',
                 paths.storage.ObjectStore(paths.Ensemble, nestable=True, has_name=True))
        # special stores
        # self.add('names', paths.storage.NameStore())

    def _initialize(self):
        # Set global attributes.
        setattr(self, 'title', 'OpenPathSampling Storage')
        setattr(self, 'ConventionVersion', '0.2')

        self.set_caching_mode('default')

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

        # TODO: Might not need to save topology

        logger.info("Saving topology")
        self.topologies.save(self.topology)

        logger.info("Create initial template snapshot")

        # Save the initial configuration
        self.snapshots.save(template)

        self.createVariable('template_idx', 'i4', 'scalar')
        self.variables['template_idx'][:] = self.snapshots.index[template]

    def _restore(self):

        self.set_caching_mode('default')

        self._restore_storages()
        self.topology = self.topologies[0]

    def set_caching_mode(self, mode='default'):
        """
        Set default values for all caches

        Parameters
        ----------
        caching : str
            One of the following values is allowed 'default', 'production',
            'analysis', 'off', 'lowmemory' and 'memtest'

        """

        available_cache_sizes = {
            'default': self.default_cache_sizes,
            'analysis': self.analysis_cache_sizes,
            'production': self.production_cache_sizes,
            'off': self.no_cache_sizes,
            'lowmemory': self.lowmemory_cache_sizes,
            'memtest': self.memtest_cache_sizes
        }

        if mode in available_cache_sizes:
            # We need cache sizes as a function. Otherwise we will reuse the same
            # caches for each storage and that will cause problems! Lots of...
            cache_sizes = available_cache_sizes[mode]()
        else:
            raise ValueError(
                "mode '" + mode + "' is not supported. Try one of " +
                str(available_cache_sizes.keys())
            )

        for store_name, caching in cache_sizes.iteritems():
            if hasattr(self, store_name):
                store = getattr(self, store_name)
                store.set_caching(caching)

    @staticmethod
    def default_cache_sizes():
        """
        Cache sizes for standard sessions for medium production and analysis.

        """

        return {
            'trajectories': WeakLRUCache(10000),
            'snapshots': WeakLRUCache(10000),
            'configurations': WeakLRUCache(10000),
            'momenta': WeakLRUCache(10000),
            'samples': WeakLRUCache(25000),
            'samplesets': False,
            'cvs': True,
            'pathmovers': True,
            'shootingpoints': WeakLRUCache(10000),
            'shootingpointselectors': True,
            'engines': True,
            'pathsimulators': True,
            'volumes': True,
            'ensembles': True,
            'pathmovechanges': False,
            'transitions': True,
            'networks': True,
            'details': False,
            'steps': WeakLRUCache(1000)
        }

    @staticmethod
    def lowmemory_cache_sizes():
        """
        Cache sizes for very low memory

        This uses even less caching than production runs. Mostly used for debugging.
        """

        return {
            'trajectories': WeakLRUCache(10),
            'snapshots': WeakLRUCache(100),
            'configurations': WeakLRUCache(10),
            'momenta': WeakLRUCache(10),
            'samples': WeakLRUCache(25),
            'samplesets': False,
            'cvs': True,
            'pathmovers': True,
            'shootingpoints': False,
            'shootingpointselectors': True,
            'engines': True,
            'pathsimulators': True,
            'volumes': True,
            'ensembles': True,
            'pathmovechanges': False,
            'transitions': True,
            'networks': True,
            'details': False,
            'steps': WeakLRUCache(10)
        }


    @staticmethod
    def memtest_cache_sizes():
        """
        Cache Sizes for memtest debugging sessions

        Memtest will cache everything weak to measure if there is some object left in
        memory that should have been disposed of.

        """
        return {
            'trajectories': WeakLRUCache(10),
            'snapshots': WeakLRUCache(10),
            'configurations': WeakLRUCache(10),
            'momenta': WeakLRUCache(10),
            'samples': WeakLRUCache(10),
            'samplesets': WeakLRUCache(10),
            'cvs': WeakLRUCache(10),
            'pathmovers': WeakLRUCache(10),
            'shootingpoints': WeakLRUCache(10),
            'shootingpointselectors': WeakLRUCache(10),
            'engines': WeakLRUCache(10),
            'pathsimulators': WeakLRUCache(10),
            'volumes': WeakLRUCache(10),
            'ensembles': WeakLRUCache(10),
            'pathmovechanges': WeakLRUCache(10),
            'transitions': WeakLRUCache(10),
            'networks': WeakLRUCache(10),
            'details': WeakLRUCache(10),
            'steps': WeakLRUCache(10)
        }

    #

    @staticmethod
    def analysis_cache_sizes():
        """
        Cache Sizes for analysis sessions

        Analysis caching is very large to allow fast processing

        """
        return {
            'trajectories': WeakLRUCache(500000),
            'snapshots': WeakLRUCache(100000),
            'configurations': WeakLRUCache(10000),
            'momenta': WeakLRUCache(1000),
            'samples': WeakLRUCache(1000000),
            'samplesets': WeakLRUCache(100000),
            'cvs': True,
            'pathmovers': True,
            'shootingpoints': WeakLRUCache(100000),
            'shootingpointselectors': True,
            'engines': True,
            'pathsimulators': True,
            'volumes': True,
            'ensembles': True,
            'pathmovechanges': WeakLRUCache(250000),
            'transitions': True,
            'networks': True,
            'details': False,
            'steps': WeakLRUCache(50000)
        }


    @staticmethod
    def production_cache_sizes():
        """
        Cache Sizes for production runs

        Production. No loading assumed, only last 1000 steps and a few other
        objects for error testing

        """
        return {
            'trajectories': WeakLRUCache(100),
            'snapshots': WeakLRUCache(100),
            'configurations': WeakLRUCache(1000),
            'momenta': WeakLRUCache(1000),
            'samples': WeakLRUCache(100),
            'samplesets': False,
            'cvs': False,
            'pathmovers': False,
            'shootingpoints': False,
            'shootingpointselectors': False,
            'engines': False,
            'pathsimulators': False,
            'volumes': False,
            'ensembles': False,
            'pathmovechanges': False,
            'transitions': False,
            'networks': False,
            'details': False,
            'steps': WeakLRUCache(10)
        }

    # No caching (so far only CVs internal storage is there)

    @staticmethod
    def no_cache_sizes():
        """
        Set cache sizes to no caching at all.

        Notes
        -----
        This is VERY SLOW and only used for debugging.
        """
        return {
            'trajectories': False,
            'snapshots': False,
            'configurations': False,
            'momenta': False,
            'samples': False,
            'samplesets': False,
            'cvs': False,
            'pathmovers': False,
            'shootingpoints': False,
            'shootingpointselectors': False,
            'engines': False,
            'pathsimulators': False,
            'volumes': False,
            'ensembles': False,
            'pathmovechanges': False,
            'transitions': False,
            'networks': False,
            'details': False,
            'steps': False
        }


class AnalysisStorage(Storage):
    """
    Open a storage in read-only and do caching useful for analysis.
    """
    def __init__(self, filename):
        """
        Parameters
        ----------
        filename : str
            The filename of the storage to be opened
        """
        super(AnalysisStorage, self).__init__(
            filename=filename,
            mode='r'
        )

        self.set_caching_mode('analysis')

        # Let's go caching
        AnalysisStorage.cache_for_analysis(self)

    @staticmethod
    def cache_for_analysis(storage):
        storage.samples.cache_all()
        storage.samplesets.cache_all()
        storage.cvs.cache_all()
        storage.volumes.cache_all()
        storage.ensembles.cache_all()
        storage.pathmovers.cache_all()
        storage.pathmovechanges.cache_all()
        storage.steps.cache_all()
#        storage.trajectories.cache_all()
