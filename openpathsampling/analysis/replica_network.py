import networkx as nx

class ReplicaNetwork(object):

    def __init__(self, repex_movers=None, ensembles=None, storage=None):
        self.analysis = { } 
        self.traces = { } 
        if repex_movers is None and ensembles is None:
            raise RuntimeError("Must define either repex_movers or ensembles")
        if ensembles is None:
            tmp_ensembles = []
            for mover in repex_movers:
                tmp_ensembles.extend(mover.ensembles)
            sort_ens = sorted(tmp_ensembles)
            ensembles = [sort_ens[0]]
            for ens in sort_ens[1:]:
                if ens != ensembles[-1]:
                    ensembles.append(ens)

        if repex_movers is None:
            repex_movers = []
            for mover in storage.pathmovers:
                if isinstance(mover, paths.ReplicaExchangeMover):
                    pass


        self.repex_movers = repex_movers
        self.ensembles = ensembles

        self.storage = storage

    def check_storage(self, storage):
        """Checks whether we have a valid storage to look at.

        If storage is given as a parameter, it overwrites the instance
        variable.
        """
        if storage != None:
            if storage != self.storage:
                self.analysis = { }
                self.traces = { } 
            self.storage = storage
        if self.storage == None:
            raise RuntimeError("No storage given for analysis")

    def analyze_exchanges(self, storage, force=False):
        self.check_storage(storage)
        if force == False and self.analysis != { }:
            return self.analysis

        pass

    def analyze_traces(self, storage, force=False):
        self.check_storage(storage)
        if force == False and self.traces != { }:
            return self.traces
        for ensemble in [s.ensemble for s in self.storage.sampleset[0]]:
            self.traces[ensemble] = trace_replicas_for_ensemble(ensemble,
                                                                self.storage)
        for replica in [s.replica for s in self.storage.sampleset[0]]:
            self.traces[replica] = trace_ensembles_for_replica(replica,
                                                               self.storage)
        return self.traces

    def diagram(self, storage=None, force=False):
        self.check_storage(storage)
        pass

    def flow(self, bottom, top, storage=None, force=False):
        self.check_storage(storage)
        traces = self.analyze_traces(storage, force)
        pass

    def one_way_trips(self, bottom, top, storage=None, force=False):
        pass

    def round_trips(self, bottom, top, storage=None, force=False):
        self.check_storage(storage)
        pass


def trace_ensembles_for_replica(replica, storage):
    trace = []
    for sset in storage.samplesets:
        trace.append(sset[replica].ensemble)
    return trace

def trace_replicas_for_ensemble(ensemble, storage):
    trace = []
    for sset in storage.samplesets:
        trace.append(sset[ensemble].replica)
    return trace

def condense_repeats(ll):
    count = 0 
    old = None
    vals = []
    for e in ll:
        if e == old:
            count += 1
        else:
            if old != None:
                vals.append((old, count))
            count = 1
            old = e
    return vals
