import openpathsampling as paths

from ..simstore.storable_functions import (
    StorableFunction, StorableFunctionConfig, wrap_numpy,
    scalarize_singletons, requires_lists_pre, requires_lists_post, Processor
)
from openpathsampling.netcdfplus import StorableNamedObject
from ..simstore.serialization_helpers import get_uuid

class CollectiveVariable(StorableFunction):
    """Wrapper around functions that map snapshots to values.

    This specializes the generic :class:`.StorableFunction` for OPS
    snapshots. In particular, it makes it so that trajectories and sampls
    are seen as lists of snapshots.
    """
    def is_scalar(self, item):
        # override is_scalar so that Trajectories and Samples are treated
        # as iterables over snapshots
        if isinstance(item, (paths.Trajectory, paths.Sample)):
            return False
        else:
            return super(CollectiveVariable, self).is_scalar(item)


class ReversibleStorableFunction(StorableFunction):
    """Wrapper around functions that don't depend on the arrow of time.

    This is separate from CoordinateFunctionCV because, in principle, both
    snapshots and trajectories can be reversed, and some functions of
    trajectories (e.g., max value of a CV) can be independent of time.
    However, it appears that trajectories do not currently implement the
    ability to store the UUID of the reversed partner.
    """
    # implementation: here we overload the _get_cached and _get_storage
    # methods to also look for the reversed snapshots
    # Current implementation (no __init__) also makes this a mix-in
    def _get_forward_and_reversed(self, func, uuid_items):
        fwd, fwd_missing = func(uuid_items)
        reversed_items = [item.reversed for item in fwd_missing.values()]
        rev_to_item = {get_uuid(item): get_uuid(item.reversed)
                       for item in reversed_items}
        rev = {get_uuid(item): item for item in reversed_items}
        bkwd, missing = func(rev)
        results = {rev_to_item[uuid]: item for uuid, item in bkwd.items()}
        missing = {rev_to_item[uuid]: item.reversed
                   for uuid, item in missing.items()}
        results.update(fwd)
        return results, missing

    def _get_cached(self, uuid_items):
        return self._get_forward_and_reversed(
            func=super(ReversibleStorableFunction, self)._get_cached,
            uuid_items=uuid_items
        )

    def _get_storage(self, uuid_items):
        return self._get_forward_and_reversed(
            func=super(ReversibleStorableFunction, self)._get_storage,
            uuid_items=uuid_items
        )

class CoordinateFunctionCV(CollectiveVariable, ReversibleStorableFunction):
    """For functions independent of arrow of time with snapshot input

    The term coordinate is used here because the most common case is that
    the CV is based on coordinates. However, this really implements a
    wrapper that can take advantage of cases where a snapshot and its
    time-reversed version always give the same function result (e.g.,
    kinetic energy would also work here).
    """
    pass  # all the implementation is in the two superclasses


class FunctionFactoryCV(CollectiveVariable):
    # TODO
    pass  # this is the GeneratorCV


class MDTrajProcessor(Processor):
    def __init__(self, topology):
        """
        Parameters
        ----------
        topology : :class:`.MDTrajTopology
            an OPS MDTraj topology wrapper
        """
        super(MDTrajProcessor, self).__init__(name='mdtraj',
                                              stage='list-pre',
                                              func=self)
        self.topology = topology

    def __call__(self, values):
        top = self.topology.mdtraj
        return paths.Trajectory(values).to_mdtraj(topology=top)

class MDTrajFunctionCV(CoordinateFunctionCV):
    def __init__(self, func, topology, func_config=None, period_min=None,
                 period_max=None, **kwargs):
        if func_config is None:
            func_config = StorableFunctionConfig([
                MDTrajProcessor(topology), wrap_numpy, scalarize_singletons,
            ])
        super(MDTrajFunctionCV, self).__init__(
            func,
            func_config=func_config,
            period_min=period_min,
            period_max=period_max,
            **kwargs
        )
        self.topology = topology
        self.mdtraj_topology = topology.mdtraj

    def to_dict(self):
        dct = super().to_dict()
        dct.update({'topology': self.topology})
        return dct


class PyEMMAFeaturizerCV(FunctionFactoryCV, CoordinateFunctionCV):
    # TODO
    pass
