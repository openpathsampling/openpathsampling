"""
Created on 06.07.2014

@author: JDC Chodera, JH Prinz
"""

import logging
import openpathsampling as paths
from openpathsampling.netcdfplus import NetCDFPlus, WeakLRUCache, ObjectStore, NamedObjectStore
import openpathsampling.engines as peng

from openpathsampling.storage import Storage

from Queue import Empty

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')



class RemoteClientObject(ObjectStore):
    def _load(self, idx):
        # No loading, only stuff that is in cache and memory
        return None

    def _save(self, obj, idx):
        # No loading, only caching and what is left in memory
        pass

    def __len__(self):
        if len(self.index) == 0:
            return 0
        else:
            return max(self.index.values())

class RemoteMasterObject(ObjectStore):

    def _load(self, idx):
        # No loading, only stuff that is in cache and memory

        result = self.storage.ask(
            "_cache_.simplifier.simplify_object(_cache_.%s[%d])" % (self.prefix, idx))

        return self.storage.simplifier.from_json(result)

    def _save(self, obj, idx):
        # No loading, only caching and what is left in memory
        pass

    def __len__(self):
        return 1000


# =============================================================================================
# OPS SPECIFIC STORAGE
# =============================================================================================

class RemoteClientStorage(Storage):
    """
    A netCDF4 wrapper to store trajectories based on snapshots of an OpenMM
    simulation. This allows effective storage of shooting trajectories
    """

    USE_FEATURE_SNAPSHOTS = True

    @property
    def n_atoms(self):
        return self.topology.n_atoms

    @property
    def n_spatial(self):
        return self.topology.n_spatial

    def __init__(self):
        """
        Create a netCDF+ storage for OPS Objects

        Parameters
        ----------
        filename : string
            filename of the netcdf file to be used or created
        mode : string, default: None
            the mode of file creation, one of ``w`` (write), ``a`` (append) or
            None, which will append any existing files.
        template : :class:`openpathsampling.Snapshot`
            a Snapshot instance that contains a reference to a Topology, the
            number of atoms and used units
        """

        self._setup_class()

        self.register_store('stores', RemoteClientObject(ObjectStore))
        self.stores.set_caching(True)

        self._create_storages()
        self._initialize()

    def _create_storages(self):
        """
        Register all Stores used in the OpenPathSampling Storage

        """

        # objects with special storages

        stores = {
            'trajectories' : paths.Trajectory,
            'snapshots' : paths.BaseSnapshot,
            'samples' : paths.Sample,
            'samplesets' : paths.SampleSet,
            'pathmovechanges' : paths.PathMoveChange,
            'steps' : paths.MCStep,
            'cvs' : paths.CollectiveVariable,
            'details' : paths.Details,
            'topologies' : peng.Topology,
            'pathmovers' : paths.PathMover,
            'shootingpointselectors' : paths.ShootingPointSelector,
            'engines' : peng.DynamicsEngine,
            'pathsimulators' : paths.PathSimulator,
            'transitions' : paths.Transition,
            'schemes' : paths.MoveScheme,
            'volumes' : paths.Volume,
            'ensembles' : paths.Ensemble
        }

        for name, obj in stores.iteritems():
            self.create_store(name, RemoteClientObject(obj))

    def _initialize(self):
        # Set global attributes.
        self.set_caching_mode('default')


class RemoteMasterStorage(Storage):
    """
    A netCDF4 wrapper to store trajectories based on snapshots of an OpenMM
    simulation. This allows effective storage of shooting trajectories
    """

    USE_FEATURE_SNAPSHOTS = True

    @property
    def n_atoms(self):
        return self.topology.n_atoms

    @property
    def n_spatial(self):
        return self.topology.n_spatial

    def __init__(self, client):
        """
        Create a netCDF+ storage for OPS Objects

        Parameters
        ----------
        filename : string
            filename of the netcdf file to be used or created
        mode : string, default: None
            the mode of file creation, one of ``w`` (write), ``a`` (append) or
            None, which will append any existing files.
        template : :class:`openpathsampling.Snapshot`
            a Snapshot instance that contains a reference to a Topology, the
            number of atoms and used units
        """

        self._setup_class()

        self.register_store('stores', RemoteMasterObject(ObjectStore))
        self.stores.set_caching(True)

        self._create_storages()
        self._initialize()

        self.client = client

        self.iopub = client.iopub_channel
        self.shell = client.shell_channel

        self.client.execute(
"""
import openpathsampling as _paths
_cache_ = _paths.storage.remote.RemoteClientStorage()
"""
        )

    def ask(self, cmd):
        uuid = self.client.execute(cmd)
        found = False
        try:
            while not found:
                msg = self.client.iopub_channel.get_msg(timeout=1.0)
                if msg['msg_type'] == 'execute_result':
                    if msg['parent_header']['msg_id'] == uuid:
                        found = True
                        result = msg['content']['data']['text/plain']

        except Empty:
            return None

        return result

    def tell(self, cmd):
        return self.client.execute(cmd)

    def get_obj(self, name):
        ss = self.ask('_cache_.simplifier.simplify(%s)' % name)
        return self.simplifier.from_json(ss)

    def put_obj(self, name, obj):

        self.tell()

    def _create_storages(self):
        """
        Register all Stores used in the OpenPathSampling Storage

        """

        # objects with special storages

        stores = {
            'trajectories' : paths.Trajectory,
            'snapshots' : paths.BaseSnapshot,
            'samples' : paths.Sample,
            'samplesets' : paths.SampleSet,
            'pathmovechanges' : paths.PathMoveChange,
            'steps' : paths.MCStep,
            'cvs' : paths.CollectiveVariable,
            'details' : paths.Details,
            'topologies' : peng.Topology,
            'pathmovers' : paths.PathMover,
            'shootingpointselectors' : paths.ShootingPointSelector,
            'engines' : peng.DynamicsEngine,
            'pathsimulators' : paths.PathSimulator,
            'transitions' : paths.Transition,
            'schemes' : paths.MoveScheme,
            'volumes' : paths.Volume,
            'ensembles' : paths.Ensemble
        }

        for name, obj in stores.iteritems():
            self.create_store(name, RemoteMasterObject(obj))

    def _initialize(self):
        # Set global attributes.
        self.set_caching_mode('default')

    def close(self):
        self.kc.stop_channels()
        self.km.shutdown_kernel()
        del self.km