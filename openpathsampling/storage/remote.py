"""
Created on 06.07.2014

@author: JDC Chodera, JH Prinz
"""

import logging
import openpathsampling as paths
from openpathsampling.netcdfplus import ObjectStore
import openpathsampling.engines as peng

from openpathsampling.storage import Storage

from Queue import Empty
import re

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
            return max(self.index.values()) + 1


class RemoteMasterObject(ObjectStore):

    def __init__(self, content_class, no_store=False):
        super(RemoteMasterObject, self).__init__(content_class)
        self.no_store = no_store

    def _load(self, idx):
        # No loading, only stuff that is in cache and memory

        result = self.storage.ask(
            "_cache_.simplifier.simplify_object(_cache_.%s[%d])" % (self.prefix, idx))

        return self.storage.simplifier.from_json(result)

    def _save(self, obj, idx):
        # No loading, only caching and what is left in memory

        if hasattr(self.storage, 'client'):
            s = self.simplifier.simplify_object(obj)
            if self.no_store:
                # this will change the store on the client-side to
                # use the appropriate reconstruction
                # on client-side we want all stores to be RemoteClientObject
                # ones. if we would use the normal push then the
                # type on the master side would be recreated
                # this should only be used for the StoreStore itself
                s['_cls'] = 'RemoteClientObject'

            self.storage.tell(
                '_str = "%s"' % s)

            self.storage.tell(
                "_obj = _cache_.simplifier.from_json(_str)")

            res = self.storage.ask(
                "_cache_.%s.save(_obj)" % self.prefix)

            return int(res)
        else:
            return None

    def save(self, obj, idx=None):
        if obj in self.index:
            # has been saved so quit and do nothing
            return self.index[obj]

        if hasattr(obj, '_idx'):
            # is a proxy of a saved object so do nothing
            return obj._idx

        if not isinstance(obj, self.content_class):
            raise ValueError(
                'This store can only store object of base type "%s". Given obj is of type "%s". You'
                'might need to use another store.' % (self.content_class, obj.__class__.__name__)
            )

        n_idx = self._save(obj, None)

        if n_idx is not None:
            self.index[obj] = n_idx

            # store the name in the cache
            if hasattr(self, 'cache'):
                self.cache[n_idx] = obj

        return n_idx

    def __len__(self):
        return 10000000


# =============================================================================================
# REMOTE SPECIFIC STORAGE
# =============================================================================================

class RemoteClientStorage(Storage):
    """
    A netCDF4 wrapper to store trajectories based on snapshots of an OpenMM
    simulation. This allows effective storage of shooting trajectories
    """

    USE_FEATURE_SNAPSHOTS = True

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


class RemoteMasterStorage(Storage):
    """
    A netCDF4 wrapper to store trajectories based on snapshots of an OpenMM
    simulation. This allows effective storage of shooting trajectories
    """

    USE_FEATURE_SNAPSHOTS = True

    ANSI_ESCAPE = re.compile(r'\x1b[^m]*m')

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

        self.client = client

        self.iopub = client.iopub_channel
        self.shell = client.shell_channel

        self.client_name = '_cache_'

        self.tell(
"""
import openpathsampling as paths
paths.netcdfplus.base.StorableObject.set_observer(True)
_cache_ = paths.storage.remote.RemoteClientStorage()
"""
        )

        self._setup_class()

        self.register_store('stores', RemoteMasterObject(ObjectStore, ))
        self.stores.set_caching(True)

        self._create_storages()
        self._initialize()

        for store, idx in self.stores.index.iteritems():
            self.tell("_cache_.register_store('%s', _cache_.stores[%d])" % (store.name, idx))
            self.tell("_cache_.stores[%d].name = '%s'" % (idx, store.name))

        self.tell("_cache_.set_caching_mode('default')")
        self.tell("_cache_.simplifier.update_class_list()")

    def _tb(self, tb):
        return '\n'.join([self.ANSI_ESCAPE.sub('', l) for l in tb])

    def ask(self, cmd):
        uuid = self.tell(cmd)
        return self.listen(uuid)

    def listen(self, uuid):
        found = False
        try:
            while not found:
                msg = self.client.iopub_channel.get_msg(timeout=1.0)
                if 'msg_id' in msg['parent_header'] and msg['parent_header']['msg_id'] == uuid:
                    if msg['msg_type'] == 'execute_result':
                        found = True
                        result = msg['content']['data']['text/plain']
                    elif msg['msg_type'] =='error':
                        found = True
                        result = None
                        error = msg['content']
                        print error['evalue']
                        print self._tb(error['traceback'])

        except Empty:
            return None

        return result

    def clear(self):
        self.iopub.get_msgs()

    def tell(self, cmd):
        return self.client.execute(cmd)

    def __getitem__(self, item):
        ss = self.ask('_cache_.simplifier.simplify(%s)' % item)
        return self.simplifier.from_json(ss)

    def __setitem__(self, key, value):
        ss = self.simplifier.simplify(value)
        self.tell('%s = _cache_.simplifier.from_json("%s")' % (key, ss))

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