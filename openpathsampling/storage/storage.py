"""
Created on 06.07.2014

@author: JDC Chodera, JH Prinz
"""

import logging

import openpathsampling as paths
from openpathsampling.netcdfplus import NetCDFPlus, WeakLRUCache, ObjectStore, \
    ImmutableDictStore, NamedObjectStore
import openpathsampling.engines as peng

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')

import time


# ==============================================================================
# OPS SPECIFIC STORAGE
# ==============================================================================

class Storage(NetCDFPlus):
    """
    A netCDF4 wrapper to store trajectories based on snapshots of an OpenMM
    simulation. This allows effective storage of shooting trajectories
    """

    @property
    def _ops_version_(self):
        version = paths.version.short_version
        return version

    USE_FEATURE_SNAPSHOTS = True

    @property
    def template(self):
        """
        Return the template snapshot from the storage

        Returns
        -------
        openpathsampling.engines.BaseSnapshot
            the initial snapshot
        """
        if self._template is None:
            self._template = self.tag['template']

        return self._template

    def clone(self, filename):
        """
        Creates a copy of the netCDF file and allows to reduce the used atoms.

        Parameters
        ----------
        filename : str
            the name of the cloned storage

        Notes
        -----
        This is mostly used to remove water but keep the data intact.

        """

        storage2 = Storage(filename=filename, template=self.template, mode='w')

        # Copy all configurations and momenta to new file in reduced form
        # use ._save instead of .save to override immutability checks etc...

        if self.reference_by_uuid:
            map(storage2.statics.save, self.statics)
            map(storage2.kinetics.save, self.kinetics)
        else:
            for obj in self.statics:
                storage2.statics._save(obj, self.statics.index[obj])
            for obj in self.kinetics:
                storage2.kinetics._save(obj, self.kinetics.index[obj])

        # All other should be copied one to one. We do this explicitly
        # although we could just copy all and exclude configurations and
        # momenta, but this seems cleaner

        for storage_name in [
            'pathmovers', 'topologies', 'networks', 'details', 'trajectories',
            'shootingpointselectors', 'engines', 'volumes',
            'samplesets', 'ensembles', 'transitions', 'steps',
            'pathmovechanges', 'samples', 'snapshots', 'pathsimulators', 'cvs'
        ]:
            self.clone_store(storage_name, storage2)

        storage2.close()

    # TODO: Need to copy cvs without caches!
    def clone_empty(self, filename):
        """
        Creates a copy of the netCDF file and replicates only the static parts

        Static parts are ensembles, volumes, engines, path movers, shooting
        point selectors.  We do not need to reconstruct collective variables
        since these need to be created again completely and then the
        necessary arrays in the file will be created automatically anyway.

        Parameters
        ----------
        filename : str
            the name of the cloned storage

        Notes
        -----
        This is mostly used to restart with a fresh file. Same setup,
        no results.
        """
        storage2 = Storage(filename=filename, template=self.template, mode='w')

        for storage_name in [
            'pathmovers', 'topologies', 'networks',
            'shootingpointselectors', 'engines', 'volumes',
            'ensembles', 'transitions', 'pathsimulators'
        ]:
            self.clone_store(storage_name, storage2)

        storage2.close()

    @property
    def n_atoms(self):
        return self.topology.n_atoms

    @property
    def n_spatial(self):
        return self.topology.n_spatial

    @property
    def topology(self):
        return self.template.topology

    def __init__(self, filename, mode=None, template=None, use_uuid=True, fallback=None):
        """
        Create a netCDF+ storage for OPS Objects

        Parameters
        ----------
        filename : string
            filename of the netcdf file to be used or created
        mode : string, default: None
            the mode of file creation, one of `'w'` (write), `'a'` (append) or
            None, which will append any existing files.
        template : :class:`openpathsampling.Snapshot`
            a Snapshot instance that contains a reference to a Topology, the
            number of atoms and used units
        """

        self._template = template
        super(Storage, self).__init__(filename, mode, use_uuid=use_uuid, fallback=fallback)

    def split_snapshots(self):
        """
        Creates two copies of the current storage. One containing trajectories and the other the rest
        """

        filename = '.'.join(self.filename.split('.')[:-1])

        filename_main = filename + '_main.nc'
        filename_data = filename + '_frames.nc'

        storage_main = Storage(filename=filename_main, template=self.template, mode='w')
        storage_data = Storage(filename=filename_data, template=self.template, mode='w')

        map(storage_data.trajectories.save, self.trajectories)
        map(storage_main.trajectories.remember, self.trajectories)

        for storage_name in [
            'steps',
            'pathmovers', 'topologies', 'networks', 'details', 'trajectories',
            'shootingpointselectors', 'engines', 'volumes',
            'samplesets', 'ensembles', 'transitions', 'pathmovechanges',
            'samples', 'pathsimulators', 'cvs'
        ]:
            map(getattr(storage_main, storage_name).save, getattr(self, storage_name))

        storage_main.close()
        storage_data.close()

    def _create_storages(self):
        """
        Register all Stores used in the OpenPathSampling Storage

        """

        # objects with special storages

        self.create_store('trajectories', paths.storage.TrajectoryStore())

        self.create_store('cvs', paths.storage.CVStore())

        self.create_store('snapshots', paths.storage.SnapshotWrapperStore())

        self.create_store('samples', paths.storage.SampleStore())
        self.create_store('samplesets', paths.storage.SampleSetStore())
        self.create_store(
            'pathmovechanges',
            paths.storage.PathMoveChangeStore()
        )
        self.create_store('steps', paths.storage.MCStepStore())


        # normal objects

        self.create_store('details', ObjectStore(paths.Details))
        self.create_store('topologies', NamedObjectStore(peng.Topology))
        self.create_store('pathmovers', NamedObjectStore(paths.PathMover))
        self.create_store('shootingpointselectors',
                          NamedObjectStore(paths.ShootingPointSelector))
        self.create_store('engines', NamedObjectStore(peng.DynamicsEngine))
        self.create_store('pathsimulators',
                          NamedObjectStore(paths.PathSimulator))
        self.create_store('transitions', NamedObjectStore(paths.Transition))
        self.create_store('networks',
                          NamedObjectStore(paths.TransitionNetwork))
        self.create_store('schemes',
                          NamedObjectStore(paths.MoveScheme))

        # stores where nestable could make sense but is disabled

        self.create_store('volumes',
                          NamedObjectStore(paths.Volume, nestable=True))
        self.create_store('ensembles',
                          NamedObjectStore(paths.Ensemble, nestable=True))

        # special stores

        self.create_store('tag', ImmutableDictStore())

    def write_meta(self):
        self.setncattr('storage_format', 'openpathsampling')
        self.setncattr('storage_version', paths.version.version)

    def _initialize(self):
        # Set global attributes.
        setattr(self, 'title', 'OpenPathSampling Storage')

        self.set_caching_mode()

        # since we want to store stuff we need to finalize stores that have not been initialized yet
        self.finalize_stores()

    def _restore(self):
        self.set_caching_mode()

    def sync_all(self):
        """
        Convenience function to use ``self.cvs`` and ``self`` at once.

        Under most circumstances, you want to sync ``self.cvs`` and ``self`` at
        the same time. This just makes it easier to do that.
        """
        self.cvs.sync_all()
        self.sync()

    def set_caching_mode(self, mode='default'):
        r"""
        Set default values for all caches

        Parameters
        ----------
        mode : str
            One of the following values is allowed `default`, `production`,
            `analysis`, `off`, `lowmemory` and `memtest`

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
            # We need cache sizes as a function. Otherwise we will reuse the
            # same caches for each storage and that will cause problems!
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

    def check_version(self):
        super(Storage, self).check_version()
        try:
            s1 = self.getncattr('storage_version')
        except AttributeError:
            logger.info(
                'Using openpathsampling Pre 1.0 version. '
                'No version detected using 0.0.0'
            )
            s1 = '0.0.0'

        s2 = self._ops_version_

        cp = self._cmp_version(s1, s2)

        if cp != 0:
            logger.info('Loading different OPS storage version. '
                        'Installed version is %s and loaded version is %s'
                        % (s2, s1))
            if cp > 0:
                logger.info('Loaded version is newer consider upgrading OPS '
                            'conda package!')
            else:
                logger.info('Loaded version is older. Should be no problem '
                            'other then missing features and information')

    @staticmethod
    def default_cache_sizes():
        """
        Cache sizes for standard sessions for medium production and analysis.

        """

        return {
            'trajectories': WeakLRUCache(10000),
            'snapshots': WeakLRUCache(10000),
            'statics': WeakLRUCache(10000),
            'kinetics': WeakLRUCache(10000),
            'samples': WeakLRUCache(25000),
            'samplesets': False,
            'cvs': True,
            'pathmovers': True,
            'shootingpointselectors': True,
            'engines': True,
            'pathsimulators': True,
            'volumes': True,
            'ensembles': True,
            'pathmovechanges': False,
            'transitions': True,
            'networks': True,
            'details': False,
            'steps': WeakLRUCache(1000),
            'topologies': True
        }

    @staticmethod
    def lowmemory_cache_sizes():
        """
        Cache sizes for very low memory

        This uses even less caching than production runs.
        Mostly used for debugging.
        """

        return {
            'trajectories': WeakLRUCache(1000),
            'snapshots': WeakLRUCache(1000),
            'statics': WeakLRUCache(10),
            'kinetics': WeakLRUCache(10),
            'samples': WeakLRUCache(25),
            'samplesets': False,
            'cvs': True,
            'pathmovers': True,
            'shootingpointselectors': True,
            'engines': True,
            'pathsimulators': True,
            'volumes': True,
            'ensembles': True,
            'pathmovechanges': False,
            'transitions': True,
            'networks': True,
            'details': False,
            'steps': WeakLRUCache(10),
            'topologies': True
        }

    @staticmethod
    def memtest_cache_sizes():
        """
        Cache Sizes for memtest debugging sessions

        Memtest will cache everything weak to measure if there is some object
        left in memory that should have been disposed of.

        """
        return {
            'trajectories': WeakLRUCache(10),
            'snapshots': WeakLRUCache(10),
            'statics': WeakLRUCache(10),
            'kinetics': WeakLRUCache(10),
            'samples': WeakLRUCache(10),
            'samplesets': WeakLRUCache(10),
            'cvs': WeakLRUCache(10),
            'pathmovers': WeakLRUCache(10),
            'shootingpointselectors': WeakLRUCache(10),
            'engines': WeakLRUCache(10),
            'pathsimulators': WeakLRUCache(10),
            'volumes': WeakLRUCache(10),
            'ensembles': WeakLRUCache(10),
            'pathmovechanges': WeakLRUCache(10),
            'transitions': WeakLRUCache(10),
            'networks': WeakLRUCache(10),
            'details': WeakLRUCache(10),
            'steps': WeakLRUCache(10),
            'topologies': WeakLRUCache(10)
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
            'statics': WeakLRUCache(10000),
            'kinetics': WeakLRUCache(1000),
            'samples': WeakLRUCache(1000000),
            'samplesets': WeakLRUCache(100000),
            'cvs': True,
            'pathmovers': True,
            'shootingpointselectors': True,
            'engines': True,
            'pathsimulators': True,
            'volumes': True,
            'ensembles': True,
            'pathmovechanges': WeakLRUCache(250000),
            'transitions': True,
            'networks': True,
            'details': False,
            'steps': WeakLRUCache(50000),
            'topologies': True
        }

    @staticmethod
    def production_cache_sizes():
        """
        Cache Sizes for production runs

        Production. No loading assumed, only last 1000 steps and a few other
        objects for error testing

        """
        return {
            'trajectories': WeakLRUCache(1000),
            'snapshots': WeakLRUCache(10000),
            'statics': WeakLRUCache(1000),
            'kinetics': WeakLRUCache(1000),
            'samples': WeakLRUCache(10000),
            'samplesets': False,
            'cvs': False,
            'pathmovers': False,
            'shootingpointselectors': False,
            'engines': False,
            'pathsimulators': False,
            'volumes': False,
            'ensembles': False,
            'pathmovechanges': False,
            'transitions': False,
            'networks': False,
            'details': False,
            'steps': WeakLRUCache(10),
            'topologies': True
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
            'statics': False,
            'kinetics': False,
            'samples': False,
            'samplesets': False,
            'cvs': False,
            'pathmovers': False,
            'shootingpointselectors': False,
            'engines': False,
            'pathsimulators': False,
            'volumes': False,
            'ensembles': False,
            'pathmovechanges': False,
            'transitions': False,
            'networks': False,
            'details': False,
            'steps': False,
            'topologies': False
        }


class AnalysisStorage(Storage):
    """
    Open a storage in read-only and do caching useful for analysis.

    """

    def __init__(self, filename):
        """
        Open a storage in read-only and do caching useful for analysis.

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
        """
        Run specific caching useful for later analysis sessions.

        Parameters
        ----------
        storage : :class:`openpathsampling.storage.Storage`
            The storage the caching should act upon.

        """

        for cv, (cv_store, cv_store_idx) in storage.snapshots.cv_list.items():
            cv_store.cache.load_max()

        stores_to_cache = ['cvs',
                           'samples',
                           'samplesets',
                           'volumes',
                           'ensembles',
                           'pathmovers',
                           'pathmovechanges',
                           'steps',
                           ]

        for store_name in stores_to_cache:
            store = getattr(storage, store_name)
            with AnalysisStorage.CacheTimer(store):
                store.cache_all()

#        storage.trajectories.cache_all()

    class CacheTimer(object):
        def __init__(self, context, store=None):
            self.store = None
            self.context = context

        def __enter__(self):
            self.time = time.time()
            return

        def __exit__(self, type, value, traceback):
            dtime = (time.time() - self.time) * 1000
            if self.store:
                logger.info('%s of store `%s` [%d] in %d ms' %
                            (self.context, self.store.name, len(self.store), dtime))
            else:
                logger.info('%s in %d ms' % (self.context, dtime))



class StorageView(object):
    """
    A View on a storage that only changes the iteration over steps.

    Can be used for bootstrapping on subsets of steps and pass this object
    to analysis routines.

    """

    class StepDelegate(object):
        """
        A delegate that will alter the ``iter()`` behaviour of the
        underlying store

        Attributes
        ----------
        store : dict-like
            the dict to be wrapped
        store : :class:`openpathsampling.netcdfplus.ObjectStore`
            a reference to an object store used

        """

        def __init__(self, store, step_range):
            self.store = store
            self.step_range = step_range

        def __iter__(self):
            for idx in self.step_range:
                yield self.store[idx]

        def __getitem__(self, item):
            return self.store[item]

        def __setitem__(self, key, value):
            self.store[key] = value

    def __init__(self, storage, step_range):
        """
        Parameters
        ----------

        storage : :class:`openpathsampling.storage.Storage`
            The storage the view is watching
        step_range : iterable
            An iterable object that species the step indices to be iterated over
            when using the view

        """
        self._storage = storage

        for store in self._storage.objects:
            setattr(self, store.prefix, store)

        self.variables = self._storage.variables
        self.units = self._storage.units
        self.vars = self._storage.vars

        self.steps = StorageView.StepDelegate(self._storage.steps, step_range)
