import openpathsampling as paths
import numpy as np
import pandas as pd
import scipy.sparse
from scipy.sparse.csgraph import reverse_cuthill_mckee
import networkx as nx

import logging
logger = logging.getLogger(__name__)

class ReplicaNetwork(object):
    """
    Analysis tool for networks of replica exchanges.
    """
    def __init__(self, repex_movers=None, ensembles=None, storage=None):
        self.analysis = { } 
        self.traces = { } 
        self.transitions = { }
        self.all_ensembles = []
        self.all_replicas = []
        self.ensemble_to_number = {}
        self.ensemble_to_string = {}
        if repex_movers is None and ensembles is None and storage is None:
            raise RuntimeError("Must define either repex_movers or ensembles")
        self.storage = storage
        if self.storage is not None:
            self.check_storage(self.storage)

        if repex_movers is None and ensembles is None:
            ensembles = self.all_ensembles

        # TODO: add support for repex_mover and ensembles
        # Currently we analyze everything in storage; this would allow us to
        # limit that analysis to a subset of moves
        #if ensembles is None:
        #    tmp_ensembles = []
        #    for mover in repex_movers:
        #        tmp_ensembles.extend(mover.ensembles)
        #    sort_ens = sorted(tmp_ensembles)
        #    ensembles = [sort_ens[0]]
        #    for ens in sort_ens[1:]:
        #        if ens != ensembles[-1]:
        #            ensembles.append(ens)

        #if repex_movers is None:
        #    repex_movers = []
        #    for mover in storage.pathmovers:
        #        if isinstance(mover, paths.ReplicaExchangeMover):
        #            pass

        #self.repex_movers = repex_movers
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
                self.ensemble_to_number = {}
                self.ensemble_to_string = {}
            self.storage = storage
        if self.storage == None:
            raise RuntimeError("No storage given for analysis")
        if self.all_replicas == [] or self.all_ensembles == []:
            reps_ens = get_all_ensembles_and_replicas(storage)
            self.all_replicas = reps_ens['replicas']
            self.all_ensembles = reps_ens['ensembles']
        if self.ensemble_to_number == {} or self.ensemble_to_string == {}:
            # set the default labels here
            self.initial_order()
            sset0 = self.storage.samplesets[0]
            labels = {e : str(sset0[e].replica) for e in self.all_ensembles}
            self.set_labels(labels)
        return self.storage

    def set_labels(self, ens2str=None):
        """
        Sets label dictionaries. Requires that you run self.initial_order
        for something first.

        Parameters
        ----------
        ens2str : dict of { Ensemble : string } pairs
            conversion of Ensemble to string label
        """
        # ensemble_to_string : returns a string value for the ensemble
        # ensemble_to_number : returns a non-neg int value (column order)
        if ens2str == None: 
            if self.ensemble_to_string == {}:
                ens2str = {k : str(self.ensemble_to_number[k]) 
                           for k in self.ensemble_to_number.keys()}
            else:
                ens2str = self.ensemble_to_string
        self.ensemble_to_string = ens2str
        self.string_to_ensemble = {self.ensemble_to_string[k] : k 
                                   for k in self.ensemble_to_string.keys()}
        self.number_to_string = {
            self.ensemble_to_number[k] : self.ensemble_to_string[k]
            for k in self.ensemble_to_number.keys()
        }
        self.string_to_number = {self.number_to_string[k] : k 
                                   for k in self.number_to_string.keys()}
        self.n_ensembles = len(self.ensemble_to_number.keys())


    def initial_order(self, index_order=None):
        """
        Sets order-based dictionaries.

        Parameters
        ----------
        index_order : list of Ensembles
            the ensembles in the desired order. Defaults order in
            self.all_ensembles
        """
        # dictionaries to be used to translate between orderings (these are
        # the defaults)
        if index_order == None:
            ensemble_to_number = {ens : self.all_ensembles.index(ens) 
                                  for ens in self.all_ensembles}
        else:
            ensemble_to_number = {ens : index_order.index(ens) 
                                  for ens in index_order}
        self.ensemble_to_number = ensemble_to_number
        self.number_to_ensemble = {ensemble_to_number[k] : k 
                                   for k in ensemble_to_number.keys()}
        self.set_labels()
        self.n_ensembles = len(self.ensemble_to_number)
        return ensemble_to_number

    def analyze_exchanges(self, storage, force=False):
        # TODO: convert this into something that yields ((repA, repB),
        # accepted): separate obtaining those tuples from adding up the
        # number of trials and acceptances -- this will make the rest of the
        # code usable for non-OPS purposes
        storage = self.check_storage(storage)
        if force == False and self.analysis != { }:
            return (self.analysis['n_trials'], self.analysis['n_accepted'])
        n_trials = 0
        self.analysis['n_trials'] = {}
        self.analysis['n_accepted'] = {}
        prev = None
        for step in storage.steps:
            pmc = step.change

            if pmc.canonical.mover is not None and pmc.canonical.mover.is_ensemble_change_mover:
                n_trials += 1
                hops = []
                for old in prev.active:
                    new = step.active
                    if old.replica != new[old.ensemble].replica:
                        # i.e., the prev and step have diff rep in same ens
                        hops.append((old.ensemble, new[old.replica].ensemble))
                for hop in hops:
                    try:
                        self.analysis['n_accepted'][hop] += 1
                    except KeyError:
                        self.analysis['n_accepted'][hop] = 1

            prev = step

        # TODO: n_trials no longer needs to be a dict, but other functions
        # expect that in output, so we return it
        for key in self.analysis['n_accepted'].keys():
            self.analysis['n_trials'][key] = n_trials
        return (self.analysis['n_trials'], self.analysis['n_accepted'])


    def analyze_traces(self, storage, force=False):
        """
        Calculates all the traces (fixed replica or fixed ensemble).

        Populates the dictionary at self.traces.
        """
        self.check_storage(storage)
        if force == False and self.traces != { }:
            return self.traces
        for ensemble in [s.ensemble for s in self.storage.steps[0].active]:
            self.traces[ensemble] = condense_repeats(
                trace_replicas_for_ensemble(ensemble, self.storage)
            )
        for replica in [s.replica for s in self.storage.steps[0].active]:
            self.traces[replica] = condense_repeats(
                trace_ensembles_for_replica(replica, self.storage)
            )
        return self.traces



    def reorder_matrix(self, matrix, index_order):
        """Return dataframe with matrix row/columns in index_order.
        
        Parameters
        ----------
        matrix : a SciPy COO sparse matrix
            input sparse matrix
        index_order : list of ensembles or None
            order to list ensembles. If None, defaults to reverse
            Cuthill-McKee order.

        Returns
        -------
        pandas.DataFrame
            dataframe with rows/columns ordered as desired
        """
        #""" matrix must be a coo_matrix (I think): do other have same `data`
        #attrib?"""
        if index_order == None:
            # reorder based on RCM from scipy.sparse.csgraph
            rcm_perm = reverse_cuthill_mckee(matrix.tocsr())
            rev_perm_dict = {k : rcm_perm.tolist().index(k) for k in rcm_perm}
            perm_i = [rev_perm_dict[ii] for ii in matrix.row]
            perm_j = [rev_perm_dict[jj] for jj in matrix.col]

            new_matrix = scipy.sparse.coo_matrix(
                (matrix.data, (perm_i, perm_j)), 
                shape=(self.n_ensembles, self.n_ensembles)
            )
            reordered_labels = [self.number_to_string[k] for k in rcm_perm]
        else:
            reordered_labels = [self.number_to_string[k] 
                                for k in self.number_to_string.keys()]
            new_matrix = matrix

        reordered = pd.DataFrame(new_matrix.todense())
        reordered.index = reordered_labels
        reordered.columns = reordered_labels
        return reordered

    def matrix_and_dataframe(self, ens_i, ens_j, data, index_order=None):
        """
        Create sparse matrix and pandas.Dataframe from ensemble data.

        Parameters
        ----------
        ens_i : list of ensembles
            the "from" ensemble
        ens_j : list of ensembles
            the "to" ensemble
        data : list of floats
            the data for the transition ensA->ensB, such that 
            matrix[ensA, ensB] = data[k] with ens_i[k]=ensA, ens_j[k]=ensB
        index_order : order of ensembles for output
            see `reorder_matrix`
        """
        self.initial_order(index_order)
        i = [self.ensemble_to_number[e] for e in ens_i]
        j = [self.ensemble_to_number[e] for e in ens_j]
        matrix = scipy.sparse.coo_matrix(
            (data, (i, j)), 
            shape=(self.n_ensembles, self.n_ensembles)
        )
        df = self.reorder_matrix(matrix, index_order)
        return (matrix, df)


    def transition_matrix(self, storage=None, index_order=None, force=False):
        """
        Create the transition matrix.

        Parameters
        ----------
        storage : paths.Storage
            input data
        index_order : list of ensembles or None
            see `reorder_matrix`
        force : bool (False)
            if True, recalculate cached values

        Returns
        -------
        pandas.DataFrame
            transition matrix
        """
        (n_try, n_acc) = self.analyze_exchanges(storage, force)
        data = []
        for k in n_try.keys():
            try:
                n_acc_k = n_acc[k]
            except KeyError:
                n_acc_k = 0
            data.append(float(n_acc_k) / n_try[k])
        ens_i, ens_j = zip(*n_try.keys())

        # this part should be the same for all matrices
        self.acceptance_matrix, df = self.matrix_and_dataframe(
            ens_i, ens_j, data, index_order
        )
        return df


    def mixing_matrix(self, storage=None, index_order=None, force=False):
        """
        Create the mixing matrix.

        Parameters
        ----------
        storage : paths.Storage
            input data
        index_order : list of ensembles or None
            see `reorder_matrix`
        force : bool (False)
            if True, recalculate cached values

        Returns
        -------
        pandas.DataFrame
            mixing matrix
        """
        (n_try, n_acc) = self.analyze_exchanges(storage, force)
        data = []
        for k in n_try.keys():
            try:
                n_acc_k = n_acc[k]
            except KeyError:
                n_acc_k = 0
            data.append(float(n_acc_k) * 0.5 / n_try[k])
        ens_ii, ens_jj = zip(*n_try.keys())
        # symmetrize
        ens_i = ens_ii + ens_jj
        ens_j = ens_jj + ens_ii
        data += data

        self.mix_matrix, df = self.matrix_and_dataframe(
            ens_i, ens_j, data, index_order
        )
        return df


    def transitions_from_traces(self, storage=None, force=False):
        """
        Calculate the transitions based on the trace of a given replica.

        This gives results normalized to *all* move types.

        Parameters
        ----------
        storage : paths.Storage
            input data
        force : bool (False)
            if True, recalculate cached values
        """
        traces = self.analyze_traces(storage, force)
        transitions = {}
        for replica in [s.replica for s in self.storage.samplesets[0]]:
            trace = traces[replica]
            hops = [(trace[i][0], trace[i+1][0]) for i in range(len(trace)-1)]

            for hop in hops:
                try:
                    transitions[hop] += 1
                except KeyError:
                    transitions[hop] = 1
        self.transitions = transitions
        return transitions


    def flow(self, bottom, top, storage=None, force=False):
        """
        Replica "flow" between ensembles `bottom` and `top`.

        Replica flow at a given ensemble measures the relative number of
        visits from replicas which has last visiting the "top" ensemble and
        those which had last visited the "bottom" ensemble. Ideal flow
        should be a straight line from 1.0 at "bottom" to 0.0 at "top".

        Parameters
        ----------
        bottom : paths.Ensemble
            "bottom" ensemble for this flow calculation
        top : paths.Ensemble
            "top" ensemble for this flow calculation
        storage : paths.Storage
            input data
        force : bool (False)
            if True, recalculate cached values


        Reference
        ---------
            Katzgraber, Trebst, Huse, and Troyer. J. Stat. Mech. 2006,
            P03018 (2006). doi:10.1088/1742-5468/2006/03/P03018
        """
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
        return {e : float(n_up[e])/n_visit[e] if n_visit[e] > 0 else 0.0
                for e in self.all_ensembles}

    def trips(self, bottom, top, storage=None, force=False):
        """
        Calculate round trips, up trips, and down trips.

        An "up" trip is the number of steps to get from ensemble `bottom` to
        ensemble `top`. A "down" trip is the reverse. A "round" trip
        consists of either an up trip followed by a down trip or vice versa.

        Parameters
        ----------
        bottom : paths.Ensemble
            ensemble to be considered the "bottom" for these trips
        top : paths.Ensemble
            ensemble to be considered the "top" for these trips
        storage : paths.Storage
            storage file
        force : bool (False)
            if True, recalculate cached

        Returns
        -------
        dict
            keys "up", "down", "round", pointing to values which are a list
            of the lengths of each trip of that type
        """
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

            rt_pairs = []
            if first_direction == 1:
                rt_pairs = zip(local_down, local_up)
            elif first_direction == -1:
                rt_pairs = zip(local_up, local_down)
            else:
                warnstr = "No first direction for replica "+str(replica)+": "
                warnstr += "Are there no 1-way trips?"
                logger.warn(warnstr)
            down_trips.extend(local_down)
            up_trips.extend(local_up)
            round_trips.extend([sum(pair) for pair in rt_pairs])

        return {'down' : down_trips, 'up' : up_trips, 'round' : round_trips}

def get_all_ensembles_and_replicas(storage, first_sampleset=True):
    """
    Retrieve all ensembles and replicas used in SampleSets

    Parameters
    ----------
    storage : paths.Storage
        storage file
    first_sampleset : bool (True)
        if True, assume that all relevant information is in the first
        SampleSet. If False, search through all saved SampleSets.

    Returns
    -------
    dict
        keys: 'ensembles', 'replicas', each containing a list
    """
    if first_sampleset:
        ensembles = [s.ensemble for s in storage.steps[0].active]
        replicas = [s.replica for s in storage.steps[0].active]
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

class ReplicaNetworkGraph(object):
    """
    Wrapper for NetworkX graph object generated by replica exchange network.

    Attributes
    ----------
    repx_network : paths.ReplicaNetwork
        replica exchange network object
    storage : paths.Storage
        file for data
    """
    def __init__(self, repx_network, storage=None):
        if storage is None:
            storage = repx_network.storage
        (n_try, n_acc) = repx_network.analyze_exchanges(storage)
        self.graph = nx.Graph()
        n_accs_adj = {}
        for k in n_try.keys():
            try:
                n_accs_adj[k] = n_acc[k]
            except KeyError:
                n_accs_adj[k] = 0

        largest_weight = max(n_accs_adj.values())
        
        for entry in n_try.keys():
            self.graph.add_edge(
                entry[0], entry[1], 
                weight=float(n_accs_adj[entry])/largest_weight
            )
        
        self.weights = [10*self.graph[u][v]['weight'] 
                        for u,v in self.graph.edges()]
        

    def draw(self, layout="graphviz"):
        """
        Lay out and draw graph.

        Parameters
        ----------
        layout : string ("graphviz")
            layout method. Default is "graphviz", which also requires
            installation of pygraphviz. 
        """
        if layout == "graphviz":
            pos = nx.graphviz_layout(self.graph)
        elif layout == "spring":
            pos = nx.spring_layout(self.graph)
        elif layout == "spectral":
            pos=nx.spectral_layout(self.graph)
        elif layout == "circular":
            pos=nx.circular_layout(self.graph)

        normal = []
        msouter = []
        minus = []
        for node in self.graph.nodes():
            if isinstance(node, paths.TISEnsemble):
                normal.append(node)
            elif isinstance(node, paths.MinusInterfaceEnsemble):
                minus.append(node)
            else:
                msouter.append(node)
        
        nx.draw_networkx_nodes(self.graph, pos, nodelist=normal, 
                               node_color='r', node_size=500)
        nx.draw_networkx_nodes(self.graph, pos, nodelist=minus, 
                               node_color='b', node_size=500)
        nx.draw_networkx_nodes(self.graph, pos, nodelist=msouter, 
                               node_color='g', node_size=500)

        nx.draw_networkx_edges(self.graph, pos, width=self.weights)


# TODO: convert these into functions that do the trace for all
# replicas/ensembles in one loop
def trace_ensembles_for_replica(replica, storage):
    """
    List of which ensemble a given replica was in at each MC step.

    Parameters
    ----------
    replica : 
        replica ID
    storage : paths.Storage
        storage file

    Returns
    -------
    list
        list of ensembles
    """
    trace = []
    storage.samples.cache_all()
    for step in storage.steps:
        sset = step.active
        trace.append(sset[replica].ensemble)
    return trace

def trace_replicas_for_ensemble(ensemble, storage):
    """
    List of which replica a given ensemble held at each MC step.

    Parameters
    ----------
    ensemble : paths.Ensemble
        selected ensemble
    storage : paths.Storage
        storage file

    Returns
    -------
    list
        list of replica IDs
    """
    trace = []
    storage.samples.cache_all()
    for step in storage.steps:
        sset = step.active
        trace.append(sset[ensemble].replica)
    return trace

def condense_repeats(ll):
    """
    Count the number of consecutive repeats in a list.

    Essentially, a way of doing `uniq -c`

    Parameters
    ----------
    ll : list
        a list

    Returns
    list of tuples
        list of 2-tuples in the format (element, repeats) where element is
        the element from the list, and repeats is the number of consecutive
        times it appeared
    """
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
    vals.append((old, count))
    return vals
