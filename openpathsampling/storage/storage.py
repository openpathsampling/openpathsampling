"""
Created on 06.07.2014

@author: JDC Chodera, JH Prinz
"""

import logging
import time

import openpathsampling as paths
from openpathsampling.netcdfplus import NetCDFPlus, WeakLRUCache, ObjectStore, \
    ImmutableDictStore, NamedObjectStore, PseudoAttributeStore

from .stores import SnapshotWrapperStore

import openpathsampling.engines as peng

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


# ==============================================================================
# OPS SPECIFIC STORAGE
# ==============================================================================

class Storage(NetCDFPlus):
    """
    Create a netCDF+ storage for OPS Objects

    A netCDF4 wrapper to store trajectories based on snapshots of an OpenMM
    simulation. This allows effective storage of shooting trajectories


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

    @property
    def _ops_version_(self):
        version = paths.version.short_version
        return version

    USE_FEATURE_SNAPSHOTS = True

    def __init__(
            self,
            filename,
            mode=None,
            template=None,
            fallback=None):

        self._template = template
        super(Storage, self).__init__(
            filename,
            mode,
            fallback=fallback)

    def _create_simplifier(self):
        super(Storage, self)._create_simplifier()
        self.simplifier.safemode = False

    def _create_storages(self):
        """
        Register all Stores used in the OpenPathSampling Storage

        """

        # objects with special storages
        self.create_store('trajectories', paths.storage.TrajectoryStore())

        # topologies might be needed fot CVs so put them here
        self.create_store('topologies', NamedObjectStore(peng.Topology))

        snapshotstore = SnapshotWrapperStore()
        self.create_store('snapshots', snapshotstore)

        self.create_store('samples', paths.storage.SampleStore())
        self.create_store('samplesets', paths.storage.SampleSetStore())
        self.create_store('movechanges',
                          paths.storage.MoveChangeStore())
        self.create_store('steps', paths.storage.MCStepStore())

        # normal objects
        self.create_store('details', ObjectStore(paths.Details))
        self.create_store('pathmovers', NamedObjectStore(paths.PathMover))
        self.create_store('shootingpointselectors',
                          NamedObjectStore(paths.ShootingPointSelector))
        self.create_store('engines', NamedObjectStore(peng.DynamicsEngine))
        self.create_store('pathsimulators',
                          paths.storage.PathSimulatorStore())
        self.create_store('transitions', NamedObjectStore(paths.Transition))
        self.create_store('networks',
                          NamedObjectStore(paths.TransitionNetwork))
        self.create_store('schemes',
                          NamedObjectStore(paths.MoveScheme))
        self.create_store('interfacesets',
                          NamedObjectStore(paths.InterfaceSet))
        self.create_store('msouters',
                          NamedObjectStore(paths.MSOuterTISInterface))

        # stores where nestable could make sense but is disabled

        self.create_store('volumes',
                          NamedObjectStore(paths.Volume, nestable=True))
        self.create_store('ensembles',
                          NamedObjectStore(paths.Ensemble, nestable=True))

        # special stores

        self.create_store('tag', ImmutableDictStore())

    @property
    def tags(self):
        return self.tag

    def write_meta(self):
        self.setncattr('storage_format', 'openpathsampling')
        self.setncattr('storage_version', paths.version.version)

    def _initialize(self):
        # Set global attributes.
        setattr(self, 'title', 'OpenPathSampling Storage')

        # backwards compatibility
        self.cvs = self.attributes

        self.set_caching_mode()

    def _restore(self):
        self.set_caching_mode()

        if hasattr(self, 'cvs'):
            logger.info('Opening an old version that handles CVs differently. '
                        'You cannot extend this file, only read it.')

            if self.mode != 'r':
                logger.info('Cannot open in append mode. Closing')
                self.close()
                raise RuntimeWarning('Closing. Cannot append incompatible '
                                     'file. You can still open readable.')
        else:
            self.cvs = self.attributes

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
            'memtest': self.memtest_cache_sizes,
            'unlimited': self.unlimited_cache_sizes
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

        for store_name, caching in cache_sizes.items():
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
            'attributes': True,
            'trajectories': WeakLRUCache(10000),
            'snapshots': WeakLRUCache(10000),
            'statics': WeakLRUCache(10000),
            'kinetics': WeakLRUCache(10000),
            'samples': WeakLRUCache(25000),
            'samplesets': WeakLRUCache(10000),
            'cvs': True,
            'pathmovers': True,
            'shootingpointselectors': True,
            'engines': True,
            'pathsimulators': True,
            'volumes': True,
            'ensembles': True,
            'movechanges': WeakLRUCache(10000),
            'transitions': True,
            'networks': True,
            'interfacesets': True,
            'schemes': True,
            'msouters': True,
            'details': WeakLRUCache(1000),
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
            'attributes': True,
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
            'movechanges': False,
            'transitions': True,
            'networks': True,
            'interfacesets': True,
            'schemes': True,
            'msouters': True,
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
            'attributes': WeakLRUCache(10),
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
            'movechanges': WeakLRUCache(10),
            'transitions': WeakLRUCache(10),
            'networks': WeakLRUCache(10),
            'interfacesets': WeakLRUCache(10),
            'schemes': WeakLRUCache(10),
            'msouters': WeakLRUCache(10),
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
            'attributes': True,
            'trajectories': WeakLRUCache(500000),
            'snapshots': WeakLRUCache(100000),
            'statics': WeakLRUCache(10000),
            'kinetics': WeakLRUCache(1000),
            'samples': WeakLRUCache(1000000),
            'samplesets': WeakLRUCache(200000),
            'cvs': True,
            'pathmovers': True,
            'shootingpointselectors': True,
            'engines': True,
            'pathsimulators': True,
            'volumes': True,
            'ensembles': True,
            'movechanges': WeakLRUCache(500000),
            'transitions': True,
            'networks': True,
            'interfacesets': True,
            'schemes': True,
            'msouters': True,
            'details': WeakLRUCache(1000),
            'steps': True,
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
            'attributes': True,
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
            'movechanges': False,
            'transitions': False,
            'networks': False,
            'interfacesets': False,
            'schemes': True,
            'msouters': False,
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
            'attributes': False,
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
            'movechanges': False,
            'transitions': False,
            'networks': False,
            'interfacesets': False,
            'schemes': False,
            'msouters': False,
            'details': False,
            'steps': False,
            'topologies': False
        }

    @staticmethod
    def unlimited_cache_sizes():
        """
        Set cache sizes to no caching at all.

        Notes
        -----
        This is VERY SLOW and only used for debugging.
        """
        return {
            'trajectories': True,
            'snapshots': True,
            'statics': True,
            'kinetics': True,
            'samples': True,
            'samplesets': True,
            'cvs': True,
            'pathmovers': True,
            'shootingpointselectors': True,
            'engines': True,
            'pathsimulators': True,
            'volumes': True,
            'ensembles': True,
            'movechanges': True,
            'transitions': True,
            'networks': True,
            'interfacesets': True,
            'schemes': True,
            'msouters': True,
            'details': True,
            'steps': True,
            'topologies': True
        }


class AnalysisStorage(Storage):
    """
    Open a storage in read-only and do caching useful for analysis.

    """

    def __init__(self, filename, caching_mode='analysis'):
        """
        Open a storage in read-only and do caching useful for analysis.

        Parameters
        ----------
        filename : str
            The filename of the storage to be opened
        caching_mode : str
            The caching mode to be used. Default is `analysis` which will
            cache lots of usually relevant object. If you have a decent
            size system and lots of memory you might want to try `unlimited`
            which will not load all objects but keep every object you load.
            This is fastest but might crash for large storages.

        """
        super(AnalysisStorage, self).__init__(
            filename=filename,
            mode='r'
        )

        self.set_caching_mode(caching_mode)

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

        with AnalysisStorage.CacheTimer('Cached all CVs'):
            for cv, cv_store in storage.snapshots.attribute_list.items():
                if cv_store:
                    cv_store.cache.load_max()

        stores_to_cache = ['cvs',
                           'trajectories',
                           'volumes',
                           'ensembles',
                           'samples',
                           'samplesets',
                           'pathmovers',
                           'movechanges',
                           'steps',
                           ]

        for store_name in stores_to_cache:
            store = getattr(storage, store_name)
            with AnalysisStorage.CacheTimer('Cache all objects', store):
                store.cache_all()

    class CacheTimer(object):
        def __init__(self, context, store=None):
            self.store = store
            self.context = context

        def __enter__(self):
            self.time = time.time()
            return

        def __exit__(self, type, value, traceback):
            dtime = (time.time() - self.time) * 1000
            if self.store:
                logger.info(
                    '%s of store `%s` [%d] in %d ms' %
                    (self.context, self.store.name, len(self.store), dtime))
            else:
                logger.info('%s in %d ms' % (self.context, dtime))
