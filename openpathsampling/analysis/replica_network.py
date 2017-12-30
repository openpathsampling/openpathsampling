import collections
import openpathsampling as paths
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
    def __init__(self, scheme, steps, replicas=None):
        if replicas is None:
            replicas = steps[0].active.replica_list()
        try:
            self.n_replicas = len(replicas)
        except TypeError:
            replicas = [replicas]
            self.n_replicas = 1
        self.replicas = replicas
        self.scheme = scheme
        self.ensembles = scheme.network.all_ensembles

        # set defaults
        self._ensemble_to_string = {}
        self.ensemble_order = scheme.network.all_ensembles

        self.traces = self._traces_from_steps(steps)
        self.transitions = self._transitions_from_traces(self.traces)
        self.analysis = self._analysis_from_steps(steps)


    def to_dict(self):
        dct = {
            'scheme': self.scheme,
            'replicas': self.replicas,
            'ensembles': self.ensembles,
            'ensemble_order': self.ensemble_order,
            'ensemble_to_string': self._ensemble_to_string,
            'traces': self.traces,
            'transitions': self.transitions,
            'analysis': self.analysis
        }
        return dct

    @classmethod
    def from_dict(cls, dct):
        obj = cls.__new__()
        obj.scheme = dct['scheme']
        obj.replicas = dct['replicas']
        obj.n_replcias = len(obj.replicas)
        obj.ensembles = dct['ensembles']
        obj.ensemble_order = dct['ensemble_order']
        obj._ensemble_to_string = dct['ensemble_to_string']
        obj.traces = dct['traces']
        obj.transitions = dct['transitions']
        obj.analysis = dct['analysis']
        return obj


    @property
    def number_to_ensemble(self):
        return {i: ens for (i, ens) in enumerate(self.ensemble_order)}

    @property
    def ensemble_to_number(self):
        return {ens: i for (i, ens) in enumerate(self.ensemble_order)}

    @property
    def string_to_ensemble(self):
        return {v: k for (k, v) in self.ensemble_to_string.items()}

    @property
    def n_ensembles(self):
        return len(self.ensemble_order)

    @property
    def number_to_string(self):
        ens2str = self.ensemble_to_string
        num2ens = self.number_to_ensemble
        return {i: ens2str[ens] for (i, ens) in num2ens.items()}

    @property
    def string_to_number(self):
        str2ens = self.string_to_ensemble
        ens2num = self.ensemble_to_number
        return {s: ens2num[ens] for (s, ens) in str2ens.items()}

    @property
    def ensemble_order(self):
        return self._ensemble_order

    @ensemble_order.setter
    def ensemble_order(self, value):
        self._ensemble_order = value
        default_names = {ens: ens.name for ens in value}
        default_names.update(self.ensemble_to_string)
        self.ensemble_to_string = default_names

    @property
    def ensemble_to_string(self):
        return self._ensemble_to_string

    @ensemble_to_string.setter
    def ensemble_to_string(self, value):
        self._ensemble_to_string.update(value)

    def _analysis_from_steps(self, steps=None):
        if steps is None:
            raise RuntimeError("No steps given to analyze!")
        n_trials = 0
        analysis = {}
        analysis['n_trials'] = {}
        analysis['n_accepted'] = {}
        prev = None
        for step in steps:
            canonical_mover = step.change.canonical.mover
            if canonical_mover and canonical_mover.is_ensemble_change_mover:
                n_trials += 1
                hops = []
                for old in prev.active:
                    new = step.active
                    if old.replica != new[old.ensemble].replica:
                        # i.e., the prev and step have diff rep in same ens
                        hops.append((old.ensemble, new[old.replica].ensemble))
                for hop in hops:
                    try:
                        analysis['n_accepted'][hop] += 1
                    except KeyError:
                        analysis['n_accepted'][hop] = 1

            prev = step

        # TODO: n_trials no longer needs to be a dict, but other functions
        # expect that in output, so we return it
        for key in analysis['n_accepted'].keys():
            analysis['n_trials'][key] = n_trials
        return analysis['n_trials'], analysis['n_accepted']


    def _traces_from_steps(self, steps):
        """
        Calculates all the traces (fixed replica or fixed ensemble).
        """
        full_traces = collections.defaultdict(list)
        for step in steps:
            for sample in step.active:
                ens = sample.ensemble
                rep = sample.replica
                full_traces[ens].append(rep)
                full_traces[rep].append(ens)

        traces = {k: condense_repeats(full_traces[k]) for k in full_traces}
        return traces


    def _transitions_from_traces(self, traces):
        """
        Calculate the transitions based on the trace of a given replica.

        This gives results normalized to *all* move types.

        Parameters
        ----------
        traces: dict
        """
        transitions = {}
        for replica in traces:
            trace = traces[replica]
            hops = [(trace[i][0], trace[i+1][0]) for i in range(len(trace)-1)]

            for hop in hops:
                try:
                    transitions[hop] += 1
                except KeyError:
                    transitions[hop] = 1
        return transitions


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
        if index_order is None:
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
        i = [self.ensemble_to_number[e] for e in ens_i]
        j = [self.ensemble_to_number[e] for e in ens_j]
        matrix = scipy.sparse.coo_matrix(
            (data, (i, j)),
            shape=(self.n_ensembles, self.n_ensembles)
        )
        df = self.reorder_matrix(matrix, index_order)
        return (matrix, df)


    def transition_matrix(self, index_order=None, force=False):
        """
        Create the transition matrix.

        Parameters
        ----------
        steps : iterable of :class:`.MCStep`
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
        (n_try, n_acc) = self.analysis
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


    def mixing_matrix(self, steps=None, index_order=None, force=False):
        """
        Create the mixing matrix.

        Parameters
        ----------
        steps : iterable of :class:`.MCStep`
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
        (n_try, n_acc) = self.analysis
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

    def flow(self, bottom, top, included_ensembles=None):
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
        steps : iterable of :class:`.MCStep`
            input data
        force : bool (False)
            if True, recalculate cached values


        References
        ----------
            Katzgraber, Trebst, Huse, and Troyer. J. Stat. Mech. 2006,
            P03018 (2006). doi:10.1088/1742-5468/2006/03/P03018
        """
        if included_ensembles is None:
            included_ensembles = self.ensembles
        traces = self.traces
        n_up = { ens : 0 for ens in self.ensembles }
        n_visit = { ens : 0 for ens in self.ensembles }
        for replica in self.replicas:
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
        as_dict =  {e : float(n_up[e])/n_visit[e] if n_visit[e] > 0 else 0.0
                    for e in included_ensembles}
        return as_dict

    def flow_pd(self, bottom, top):
        flow_dict = self.flow(bottom, top, steps, force)
        inverted = {flow_dict[k]: k for k in flow_dict.keys()}
        sorted_inverted = list(reversed(sorted(inverted.keys())))
        re_keyed = {k: sorted_inverted[k]
                    for k in range(len(sorted_inverted))}
        return pd.Series(re_keyed)


    def trips(self, bottom, top):
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
        steps : iterable of :class:`.MCStep`
            input data
        force : bool (False)
            if True, recalculate cached

        Returns
        -------
        dict
            keys "up", "down", "round", pointing to values which are a list
            of the lengths of each trip of that type
        """
        traces = self.traces
        down_trips = []
        up_trips = []
        round_trips = []
        for replica in self.replicas:
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


class ReplicaNetworkGraph(object):
    """
    Wrapper for NetworkX graph object generated by replica exchange network.

    Attributes
    ----------
    repx_network : paths.ReplicaNetwork
        replica exchange network object
    """
    def __init__(self, repx_network):
        (n_try, n_acc) = repx_network.analysis
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
def trace_ensembles_for_replica(replica, steps):
    """
    List of which ensemble a given replica was in at each MC step.

    Parameters
    ----------
    replica :
        replica ID
    steps : iterable of :class:`.MCStep`
        input data

    Returns
    -------
    list
        list of ensembles
    """
    return [s.active[replica].ensemble for s in steps]


def trace_replicas_for_ensemble(ensemble, steps):
    """
    List of which replica a given ensemble held at each MC step.

    Parameters
    ----------
    ensemble : paths.Ensemble
        selected ensemble
    steps : iterable of :class:`.MCStep`
        MC steps to analyze (nonsense if not contiguous

    Returns
    -------
    list
        list of replica IDs
    """
    trace = []
    for step in steps:
        sset = step.active
        trace.append(sset[ensemble].replica)
    return trace


def condense_repeats(ll, use_is=True):
    """
    Count the number of consecutive repeats in a list.

    Essentially, a way of doing `uniq -c`

    Parameters
    ----------
    ll : list
        a list

    Returns
    -------
    list of tuples
        list of 2-tuples in the format (element, repeats) where element is
        the element from the list, and repeats is the number of consecutive
        times it appeared
    """
    count = 0
    old = None
    vals = []
    for e in ll:
        if (use_is and e is old) or (not use_is and e == old):
            count += 1
        else:
            if old is not None:
                vals.append((old, count))
            count = 1
            old = e
    vals.append((old, count))
    return vals
