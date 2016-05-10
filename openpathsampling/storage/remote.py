"""
Created on 06.07.2014

@author: JDC Chodera, JH Prinz
"""

import logging
import openpathsampling as paths
from openpathsampling.netcdfplus import ObjectStore
import openpathsampling.engines as peng

from openpathsampling.storage import Storage
from openpathsampling.netcdfplus.dictify import UUIDObjectJSON

from Queue import Empty
import re

from uuid import UUID

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

    def _set_uuid(self, idx, uuid):
        pass

    @property
    def uuid_idx(self):
        return self._uuid_idx


class RemoteMasterObject(ObjectStore):

    length = 0L

    def __init__(self, content_class, no_store=False):
        super(RemoteMasterObject, self).__init__(content_class)
        self.no_store = no_store

    def __len__(self):
        return self.length

    @property
    def uuid_idx(self):
        return self._uuid_idx

    def load(self, idx):
        """
        Returns an object from the storage.

        Parameters
        ----------
        idx : int
            the integer index of the object to be loaded

        Returns
        -------
        :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the loaded object
        """

        n_idx = -1

        if self.reference_by_uuid and type(idx) is UUID:
            # we want to load by uuid and it was not in cache.
            if str(idx) in self.uuid_idx:
                n_idx = int(self.uuid_idx[str(idx)])

                # if it is in the cache, return it
                try:
                    obj = self.cache[n_idx]
                    logger.debug('Found IDX #' + str(idx) + ' in cache. Not loading!')
                    return obj

                except KeyError:
                    pass

        else:
            raise ValueError(
                'indices of type "%s" are not allowed in named storage (only UUID)' % type(idx).__name__
            )


        logger.debug('Calling load object of type ' + self.content_class.__name__ + ' and IDX #' + str(idx))

        # this loader uses UUID not an integer index to load from remote kernel
        obj = self._load(idx)

        # let's add the loaded object to this cache also and for that we need an index

        n_idx = self.free()
        self.index[obj] = n_idx

        if obj is not None:
            # update cache there might have been a change due to naming
            self.cache[n_idx] = obj

            if self.reference_by_uuid:
                setattr(obj, '__uuid__', idx)
                # make sure that you cannot change the uuid of loaded objects

                # finally store the uuid of a uuidd object in cache
                self._update_uuid_in_cache(obj.__uuid__, n_idx)

        return obj

    def _load(self, idx):
        # No loading, only stuff that is in cache and memory

        uuid = str(idx)

        result = self.storage.ask(
            "_cache_.simplifier.simplify_object(_cache_.%s[paths.storage.remote.UUID('%s')])" % (self.prefix, uuid))

        # check_uuid = self.storage.ask(
        #     "_cache_.%s[paths.storage.remote.UUID('%s')].__uuid__" % (self.prefix, uuid)
        # )
        #
        # assert(uuid == check_uuid[6:-2])

        obj = self.storage.simplifier.from_json(result)
        obj.__uuid__ = UUID(uuid)

        return obj

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

            self.storage.tell(
                "_obj.__uuid__ = paths.storage.remote.UUID('%s')" % obj.__uuid__
            )

            _ = self.storage.ask(
                "_cache_.%s.save(_obj)" % self.prefix)

        self.length = max(self.length, idx + 1)

    def _set_uuid(self, idx, uuid):
        pass


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

        self.reference_by_uuid = True
        self.simplifier = UUIDObjectJSON(self)

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

        self.reference_by_uuid = True

        self.tell(
"""
import openpathsampling as paths
# paths.netcdfplus.base.StorableObject.set_observer(True)
_cache_ = paths.storage.remote.RemoteClientStorage()
"""
        )

        self.simplifier = UUIDObjectJSON(self)

        self._setup_class()

        self.register_store('stores', RemoteMasterObject(ObjectStore, True))
        self.stores.set_caching(True)

        self._create_storages()
        self._initialize()

        for uuid, idx in self.stores._uuid_idx.iteritems():
            store = self.stores.cache[idx]
            self.tell("_cache_.register_store('%s', _cache_.stores[paths.storage.remote.UUID('%s')])" % (store.name, uuid))
            self.tell("_cache_.%s.name = '%s'" % (store.name, store.name))

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
#                print msg
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