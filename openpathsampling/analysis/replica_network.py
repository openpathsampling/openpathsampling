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
        self.all_ensembles = []
        self.all_replicas = []

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
                self.all_replicas = []
                self.all_ensembles = []
            self.storage = storage
        if self.storage == None:
            raise RuntimeError("No storage given for analysis")
        if self.all_replicas == [] or self.all_ensembles == []:
            reps_ens = get_all_ensembles_and_replicas(storage)
            self.all_replicas = reps_ens['replicas']
            self.all_ensembles = reps_ens['ensembles']


    def analyze_exchanges(self, storage, force=False):
        self.check_storage(storage)
        if force == False and self.analysis != { }:
            return self.analysis
        # TODO: this generates two matrices: naccepted and ntrials. Each
        # should be represented as a sparse upper triangular matrix.
        pass

    def analyze_traces(self, storage, force=False):
        self.check_storage(storage)
        if force == False and self.traces != { }:
            return self.traces
        for ensemble in [s.ensemble for s in self.storage.sampleset[0]]:
            self.traces[ensemble] = condense_repeats(
                trace_replicas_for_ensemble(ensemble, self.storage)
            )
        for replica in [s.replica for s in self.storage.sampleset[0]]:
            self.traces[replica] = condense_repeats(
                trace_ensembles_for_replica(replica, self.storage)
            )
        return self.traces

    def transition_matrix(self, storage=None, index_order=None, force=False):
        (nacc, ntry) = self.analysis_exchanges(storage, force)
        # TODO: convert it to a pandas dataframe and return it


    def diagram(self, storage=None, force=False):
        (nacc, ntry) = self.analysis_exchanges(storage, force)
        # TODO: make this into a networkx diagram. It would be really nice
        # if a given interface set could be forced to be collinear

    def flow(self, bottom, top, storage=None, force=False):
        traces = self.analyze_traces(storage, force)
        pass

    def one_way_trips(self, bottom, top, storage=None, force=False):
        traces = self.analyze_traces(storage, force)
        pass

    def round_trips(self, bottom, top, storage=None, force=False):
        traces = self.analyze_traces(storage, force)
        pass

def get_all_ensembles_and_replicas(storage, first_sampleset=True):
    if first_sampleset:
        ensembles = [s.ensemble for s in storage.sampleset[0]]
        replicas = [s.replica for s in storage.sampleset[0]]
    else:
        # This approach uses dicts so we don't have to hunt for the key; the
        # value assigned is arbitrarily 1. Still has to loop over
        # nsets*nsamples, but that's better than nsets*nsamples*nensembles
        ensembles_dict = {}
        replicas_dict = {}
        for sset in storage.sampleset:
            for s in sset:
                ensembles[s.ensemble] = 1
                replicas[s.replica] = 1
        ensembles = ensembles_dict.keys()
        replicas = replicas_dict.keys()
    return { 'ensembles' : ensembles, 'replicas' : replicas }


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
