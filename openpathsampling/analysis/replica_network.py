import openpathsampling as paths
import numpy as np
import pandas as pd
import scipy.sparse
from scipy.sparse.csgraph import reverse_cuthill_mckee
import networkx as nx



class ReplicaNetwork(object):

    def __init__(self, repex_movers=None, ensembles=None, storage=None):
        self.analysis = { } 
        self.traces = { } 
        self.all_ensembles = []
        self.all_replicas = []
        if repex_movers is None and ensembles is None and storage is None:
            raise RuntimeError("Must define either repex_movers or ensembles")
        self.storage = storage
        if self.storage is not None:
            self.check_storage(self.storage)
        if repex_movers is None and ensembles is None:
            ensembles = self.all_ensembles

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
        # TODO: convert this into something that yields ((repA, repB),
        # accepted): separate obtaining those tuples from adding up the
        # number of trials and acceptances -- this will make the rest of the
        # code usable for non-OPS purposes
        self.check_storage(storage)
        if force == False and self.analysis != { }:
            return (self.analysis['n_trials'], self.analysis['n_accepted'])
        self.analysis['n_trials'] = {}
        self.analysis['n_accepted'] = {}
        for pmc in storage.pathmovechanges:
            for delta in pmc:
                if isinstance(delta.mover, paths.ReplicaExchangeMover):
                    if len(delta.trials) == 2:
                        ens1 = delta.trials[0].ensemble
                        ens2 = delta.trials[1].ensemble
                    else:
                        print "RepEx mover with n_trials != 2"
                        try:
                            # TODO: this hack for minus should not be
                            # necessary; although we may have to hack minus
                            # to be cleaner
                            ens1 = delta.mover.innermost_ensemble
                            ens2 = delta.mover.minus_ensemble
                        except:
                            raise RuntimeWarning("RepEx mover with n_trials != 2")
                    try:
                        self.analysis['n_trials'][(ens1, ens2)] += 1
                    except KeyError:
                        self.analysis['n_trials'][(ens1, ens2)] = 1
                    if delta.accepted:
                        try:
                            self.analysis['n_accepted'][(ens1, ens2)] += 1
                        except KeyError:
                            self.analysis['n_accepted'][(ens1, ens2)] = 1

        return (self.analysis['n_trials'], self.analysis['n_accepted'])

    def analyze_traces(self, storage, force=False):
        self.check_storage(storage)
        if force == False and self.traces != { }:
            return self.traces
        for ensemble in [s.ensemble for s in self.storage.samplesets[0]]:
            self.traces[ensemble] = condense_repeats(
                trace_replicas_for_ensemble(ensemble, self.storage)
            )
        for replica in [s.replica for s in self.storage.samplesets[0]]:
            self.traces[replica] = condense_repeats(
                trace_ensembles_for_replica(replica, self.storage)
            )
        return self.traces

    def initial_order(self, index_order):
        # dictionaries to be used to translate between orderings (these are
        # the defaults)
        if index_order == None:
            ensemble_to_number = {ens : self.all_ensembles.index(ens) 
                                  for ens in self.all_ensembles}
        else:
            ensemble_to_number = {ens : index_order.index(ens) 
                                  for ens in index_order}
        return ensemble_to_number


    def transition_matrix(self, storage=None, index_order=None, force=False):
        (n_try, n_acc) = self.analyze_exchanges(storage, force)
        ensemble_to_number = self.initial_order(index_order)
        number_to_ensemble = {ensemble_to_number[k] : k for 
                              k in ensemble_to_number.keys()}
        n_ensembles = len(ensemble_to_number)
        data = [float(n_acc[k]) / n_try[k] for k in n_try.keys()]
        ens_i, ens_j = zip(*n_try.keys())
        i = [ensemble_to_number[e] for e in ens_i]
        j = [ensemble_to_number[e] for e in ens_j]
        acc_matrix = scipy.sparse.coo_matrix(
            (data, (i, j)), 
            shape=(n_ensembles, n_ensembles)
        )
        # TODO clean these up: maybe move labels to elsewhere?
        sset0 = storage.samplesets[0]
        labels = {k : sset0[number_to_ensemble[k]].replica 
                  for k in number_to_ensemble.keys()}

        df = self.reorder_matrix(acc_matrix, labels, index_order)
        return df
        # TODO: convert it to a pandas dataframe and return it


    def reorder_matrix(self, matrix, number_to_label, index_order):
        """ matrix must be a coo_matrix (I think): do other have same `data`
        attrib?"""
        n_ensembles = len(number_to_label)
        if index_order == None:
            # reorder based on RCM from scipy.sparse.csgraph
            rcm_perm = reverse_cuthill_mckee(matrix.tocsr())
            perm_i = [rcm_perm[ii] for ii in matrix.row]
            perm_j = [rcm_perm[jj] for jj in matrix.col]

            new_matrix = scipy.sparse.coo_matrix(
                (matrix.data, (perm_i, perm_j)), 
                shape=(n_ensembles, n_ensembles)
            )
            reordered_labels = [number_to_label[k] for k in rcm_perm]
        else:
            reordered_labels = [number_to_label[k] 
                                for k in number_to_label.keys()]
            new_matrix = acc_matrix

        reordered = pd.DataFrame(new_matrix.todense())
        reordered.index = reordered_labels
        reordered.columns = reordered_labels
        return reordered



    def mixing_matrix(self, storage=None, index_order=None, force=False):
        (n_try, n_acc) = self.analyze_exchanges(storage, force)
        # TODO: if the transition matrix works, this just involves
        # modifying the input i, j, data

    def diagram(self, storage=None, force=False):
        (nacc, ntry) = self.analyze_exchanges(storage, force)
        # TODO: make this into a networkx diagram. It would be really nice
        # if a given interface set could be forced to be collinear

    def flow(self, bottom, top, storage=None, force=False):
        traces = self.analyze_traces(storage, force)
        n_up = { ens : 0 for ens in self.all_ensembles }
        n_visit = { ens : 0 for ens in self.all_ensembles } 
        for replica in self.all_replicas:
            trace = self.traces[replica]
            direction = 0
            for (loc, count) in trace:
                if loc == top:
                    direction = -1
                elif loc == bottom:
                    direction = +1
                if direction != 0:
                    n_visit[loc] += count
                if direction == 1:
                    n_up[loc] += count
        self._flow_up = n_up
        self._flow_count = n_visit
        return {e : float(n_up[e])/n_visit[e] for e in self.all_ensembles}

    def trips(self, bottom, top, storage=None, force=False):
        traces = self.analyze_traces(storage, force)
        down_trips = []
        up_trips = []
        round_trips = []
        for replica in self.all_replicas:
            trace = traces[replica]
            direction = None
            trip_counter = 0
            first_direction = None
            local_down = []
            local_up = []
            for (loc, count) in trace:
                if loc == top and direction != +1:
                    direction = +1
                    if trip_counter > 0:
                        local_up.append(trip_counter)
                    trip_counter = 0
                elif loc == bottom and direction != -1:
                    direction = -1
                    if trip_counter > 0:
                        local_down.append(trip_counter)
                    trip_counter = 0
                if direction is not None:
                    if first_direction is None:
                        first_direction = direction
                    trip_counter += count

            if first_direction == 1:
                rt_pairs = zip(local_down, local_up)
            elif first_direction == -1:
                rt_pairs = zip(local_up, local_down)
            else:
                raise RuntimeWarning(
                    "No first direction identified: Are there no 1-way trips?"
                )
            down_trips.extend(local_down)
            up_trips.extend(local_up)
            round_trips.extend([sum(pair) for pair in rt_pairs])

        return {'down' : down_trips, 'up' : up_trips, 'round' : round_trips}

def get_all_ensembles_and_replicas(storage, first_sampleset=True):
    if first_sampleset:
        ensembles = [s.ensemble for s in storage.samplesets[0]]
        replicas = [s.replica for s in storage.samplesets[0]]
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
    storage.samples.cache_all()
    for sset in storage.samplesets:
        trace.append(sset[replica].ensemble)
    return trace

def trace_replicas_for_ensemble(ensemble, storage):
    trace = []
    storage.samples.cache_all()
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
