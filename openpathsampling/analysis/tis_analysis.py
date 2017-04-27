import collections
import openpathsampling as paths
from openpathsampling.netcdfplus import StorableNamedObject
from openpathsampling.numerics import LookupFunction
import pandas as pd

def steps_to_weighted_trajectories(steps, ensembles):
    """Bare function to convert to the weighted trajs dictionary.

    This prepares data for the faster analysis format. This preparation only
    need to be done once, and it will cover a lot of the analysis cases.
    """
    results = {e: collections.Counter() for e in ensembles}

    # loop over blocks # TODO: add blocksize parameter, test various sizes
    block_steps = steps
    block = collections.defaultdict(list)
    for step in block_steps:
        for ens in ensembles:
            block[ens].append(step.active[ens].trajectory)

    block_counter = {e: collections.Counter(block[e]) for e in ensembles}

    for e in results:
        results[e] += block_counter[e]

    return results


class MultiEnsembleSamplingAnalyzer(StorableNamedObject):
    """Abstract class for getting statistics for MC steps sampling multiple
    ensembles."""
    def calculate(self, steps, ensembles):
        weighted_trajs = steps_to_weighted_trajectories(steps, ensembles)
        return self.from_weighted_trajectories(weighted_trajs)

    def from_weighted_trajectories(self, input_dict):
        raise NotImplementedError

    def combine_results(self, result_1, result_2):
        # to be used to simplify parallelization
        raise NotImplementedError

######## CALCULATING THE FLUX

#class MinusMoveFlux(MultiEnsembleSamplingAnalyzer):
    #def from_weighted_trajectories(self, input_dict):
        #pass

class DictFlux(MultiEnsembleSamplingAnalyzer):
    """Pre-calculated flux, provided as a dict.
    """
    def __init__(self, flux_dict):
        super(DictFlux, self).__init__()
        self.flux_dict = flux_dict

    def calculate(self, steps):
        return self.flux_dict

    def from_weighted_trajectories(self, input_dict):
        return self.flux_dict

    def combine_results(self, result_1, result_2):
        if result_1 != result_2:
            raise RuntimeError("Combining results from different DictFlux")
        return self

########## GENERAL HISTOGRAMMING

class EnsembleHistogrammer(MultiEnsembleSamplingAnalyzer):
    """
    Generic code to calculate the properly weighted histograms of trajectory
    properties per ensemble.
    """
    def __init__(self, ensembles, f, hist_parameters):
        self.ensembles = ensembles
        self.f = f
        self.hist_parameters = hist_parameters
        self.hists = {e: paths.numerics.Histogram(**self.hist_parameters)
                      for e in self.ensembles}

    def from_weighted_trajectories(self, input_dict):
        for ens in self.hists:
            trajs = input_dict[ens].keys()
            weights = input_dict[ens].values()
            data = [self.f(traj) for traj in trajs]
            self.hists[ens].histogram(data, weights)
        return self.hists


class PathLengthHistogrammer(EnsembleHistogrammer):
    """Histogramming path length distribution"""
    def __init__(self, ensembles, hist_parameters=None):
        if hist_parameters is None:
            pass  # set defaults
        super(PathLengthHistogrammer, self).__init__(
            ensembles=ensembles,
            f=lambda t: len(t),
            hist_parameters=hist_parameters
        )

############### HISTOGRAMMING MAX LAMBDA

class FullHistogramMaxLambdas(EnsembleHistogrammer):
    """Histogramming the full max-lambda function (one way of getting TCP)

    One of these per transition.
    """
    def __init__(self, transition, hist_parameters, max_lambda_func=None):
        self.transition = transition
        if max_lambda_func is None:
            max_lambda_func = lambda t: max(transition.interfaces.cv(t))
            #max_lambda_func = transition.interfaces.max_cv  # TODO traj-cv
        self.lambdas = {e: l for (e, l) in zip(transition.ensembles,
                                               transition.interfaces.lambdas)}
        super(FullHistogramMaxLambdas, self).__init__(
            ensembles=transition.ensembles,
            f=max_lambda_func,
            hist_parameters=hist_parameters
        )

#class PerEnsembleMaxLambdas(EnsembleHistogrammer):
    # TODO: this just maps the count to the ensemble, not the full histogram
    #def __init__(self, transition):
        #interfaces_lambdas = transition.interfaces.lambdas

class TotalCrossingProbability(MultiEnsembleSamplingAnalyzer):
    def __init__(self, max_lambda_calc, combiner=None):
        self.max_lambda_calc = max_lambda_calc
        self.transition = max_lambda_calc.transition
        if combiner is None:
            lambdas = self.transition.interfaces.lambdas
            combiner = paths.numerics.WHAM(interfaces=lambdas)
        self.combiner = combiner

    def from_weighted_trajectories(self, input_dict):
        hists = self.max_lambda_calc.from_weighted_trajectories(input_dict)
        return self.from_ensemble_histograms(hists)

    def from_ensemble_histograms(self, hists):
        tcp_results = {}
        input_hists = [hists[ens] for ens in self.transition.ensembles]
        df = paths.numerics.histograms_to_pandas_dataframe(
            input_hists,
            fcn="reverse_cumulative"
        ).sort_index(axis=1)
        tcp = self.combiner.wham_bam_histogram(df).to_dict()
        return LookupFunction(tcp.keys(), tcp.values())


class ConditionalTransitionProbability(MultiEnsembleSamplingAnalyzer):
    def __init__(self, ensembles, states):
        self.states = states
        self.ensembles = ensembles

    def from_weighted_trajectories(self, input_dict):
        ctp = {}
        for ens in self.ensembles:
            acc = collections.Counter()
            n_try = sum(input_dict[ens].values())
            final_frames = [traj.get_as_proxy(-1) for
                            traj in input_dict[ens].keys()]
            weights = input_dict[ens].values()
            for (f, w) in zip(final_frames, weights):
                local = collections.Counter({s: w for s in self.states
                                             if s(f)})
                acc += local

            ctp[ens] = {s : float(acc[s]) / n_try for s in acc.keys()}
            # TODO: add logging to report here
        return ctp

class StandardTransitionProbability(MultiEnsembleSamplingAnalyzer):
    def __init__(self, transition, tcp_method, ctp_method):
        self.transition = transition
        self.tcp_method = tcp_method
        self.ctp_method = ctp_method
        self.final_state = self.transition.stateB
        self.outermost_ensemble = self.transition.ensembles[-1]
        interfaces = transition.interfaces
        self.outermost_lambda = interfaces.get_lambda(interfaces[-1])

    def from_weighted_trajectories(self, input_dict):
        tcp = self.tcp_method.from_weighted_trajectories(input_dict)
        ctp = self.ctp_method.from_weighted_trajectories(input_dict)
        return self.from_intermediate_results(tcp, ctp)

    def from_intermediate_results(self, tcp, ctp):
        outermost_ctp = ctp[self.outermost_ensemble][self.final_state]
        tcp_at_outermost = tcp(self.outermost_lambda)
        # TODO: log things here
        return outermost_ctp * tcp_at_outermost


class TransitionDictResults(StorableNamedObject):
    # allows you to use analysis transition, 2-tuple of states, or sampling
    # transition as the key to retrieve the stored results
    def __init__(self, results_dict, network):
        self.results_dict = results_dict
        self.network = network

    def __getitem__(self, key):
        if key in self.network.sampling_transitions:
            key = self.network.sampling_to_analysis(key)
        try:
            result = self.results_dict[key]
        except KeyError:
            result = self.results_dict[self.network.transitions[key]]
        return result

    def to_pandas(self, order=None):
        key_map = lambda key: key.name
        keys = self.results_dict.keys()
        idx_vols = [k[0] for k in keys]
        col_vols = [k[1] for k in keys]
        if order is None:
            order = set(idx_vols + col_vols)
        index = [key_map(k) for k in order if k in idx_vols]
        columns = [key_map(k) for k in order if k in col_vols]
        result = pd.DataFrame(index=index, columns=columns)
        for k in keys:
            result.set_value(key_map(k[0]), key_map(k[1]),
                             self.results_dict[k])
        return result

    def __str__(self):
        return self.to_pandas().__str__()

    def __repr__(self):
        return self.to_pandas().__repr__()


class TISAnalysis(StorableNamedObject):
    """
    Generic class for TIS analysis. One of these for each network.

    Parameters
    ----------
    network : :class:`.TransitionNetwork`
    flux_method : flux calculation method
    transition_probability_methods : dict of :class:`.Transition` to method
    """
    def __init__(self, network, flux_method, transition_probability_methods):
        self.network = network
        self.transitions = network.transitions
        self.flux_method = flux_method
        self.transition_probability_methods = transition_probability_methods
        self.results = {}

    def calculate(self, steps):
        self.results = {}
        weighted_trajs = steps_to_weighted_trajectories(
            steps,
            self.network.sampling_ensembles
        )
        self.from_weighted_trajectories(weighted_trajs)

    def from_weighted_trajectories(self, input_dict):
        flux_m = self.flux_method
        tp_m = self.transition_probability_methods
        fluxes = flux_m.from_weighted_trajectories(input_dict)
        # dict of transition to transition probability
        trans_prob = {t: tp_m[t].from_weighted_trajectories(input_dict)
                      for t in tp_m.keys()}
        self.results['flux'] = fluxes
        self.results['transition_probability'] = TransitionDictResults(
            {(t.stateA, t.stateB) : trans_prob[t] for t in trans_prob},
            self.network
        )

        rates = {}
        for (trans, transition_probability) in trans_prob.iteritems():
            trans_flux = fluxes[(trans.stateA, trans.interfaces[0])]
            rates[(trans.stateA, trans.stateB)] = \
                    trans_flux * transition_probability

        self.results['rate'] = TransitionDictResults(rates, self.network)
        return self.results

    def _access_cached_result(self, key):
        try:
            return self.results[key]
        except KeyError:
            raise AttributeError("Can't access results for '" + key
                                 + "' until analysis is performed")

    @property
    def flux_matrix(self):
        return self._access_cached_result('flux')

    def flux(self, from_state, through_interface=None):
        fluxes = self._access_cached_result('flux')
        if through_interface is None:
            through_interface = from_state

        return fluxes[(from_state, through_interface)]

    def state_fluxes(self, from_state):
        fluxes = self._access_cached_result('flux')
        state_keys = [k for k in fluxes.keys() if k[0] == from_state]
        return {k: fluxes[k] for k in state_keys}

    @property
    def transition_probability_matrix(self):
        return self._access_cached_result('transition_probability')

    def transition_probability(self, from_state, to_state):
        trans_probs = self._access_cached_result('transition_probability')
        return trans_probs[(from_state, to_state)]

    def rate_matrix(self, steps=None):
        if steps is not None:
            self.calculate(steps)
        return self._access_cached_result('rate')

    def rate(self, from_state, to_state):
        return self._access_cached_result('rate')[(from_state, to_state)]


class StandardTISAnalysis(TISAnalysis):
    def __init__(self, network, steps=None, flux_method=None,
                 ctp_method=None, max_lambda_calcs=None, combiners=None):
        # NOTE: each of flux, ctp, tcp refer to the methods used; in
        # principle, these should have the option of being provided as a
        # single example (to be applied to all) or as a dict showing which
        # to apply to which analysis transition

        # set default analysis behaviors
        if flux_method is None:
            flux_pairs = None  # TODO
            flux_method = MinusMoveFlux(flux_pairs)

        if max_lambda_calcs is None:
            raise RuntimeError("Must set either max_lambda_calcs "
                               " in StandardTISAnalysis")
        max_lambda_calc_dict = {}
        for (transition, calc) in max_lambda_calcs.iteritems():
            if isinstance(calc, EnsembleHistogrammer):
                max_lambda_calc_dict[transition] = calc
            elif isinstance(calc, dict):
                max_lambda_calc_dict[transition] = FullHistogramMaxLambdas(
                    transition=transition,
                    hist_parameters=calc
                )
        if combiners is None:
            combiners = {
                transition.interfaces:
                paths.numerics.WHAM(interfaces=transition.interfaces.lambdas)
                for transition in network.sampling_transitions
            }
        self.tcp_methods = {
            transition: TotalCrossingProbability(
                max_lambda_calc=max_lambda_calc_dict[transition],
                combiner=combiners[transition.interfaces]
            )
            for transition in network.sampling_transitions
        }
        if ctp_method is None:
            outermost_ensembles = [t.ensembles[-1]
                                   for t in network.sampling_transitions]
            self.ctp_method = \
                    ConditionalTransitionProbability(
                        ensembles=outermost_ensembles,
                        states=network.all_states
                    )

        trans_prob_methods = {
            trans: StandardTransitionProbability(
                transition=trans,
                tcp_method=self.tcp_methods[trans],
                ctp_method=self.ctp_method
            )
            for trans in network.sampling_transitions
        }

        super(StandardTISAnalysis, self).__init__(
            network=network,
            flux_method=flux_method,
            transition_probability_methods=trans_prob_methods
        )

        if steps is not None:
            self.calculate(steps)

    def from_weighted_trajectories(self, input_dict):
        # calculate fluxes
        flux_m = self.flux_method
        fluxes = flux_m.from_weighted_trajectories(input_dict)
        self.results['flux'] = fluxes

        # calculate the max_lambda hists
        max_lambda_calcs = [tcp_m.max_lambda_calc
                            for tcp_m in self.tcp_methods.values()]
        max_lambda_hists = {}
        for calc in max_lambda_calcs:
            max_lambda_hists.update(
                calc.from_weighted_trajectories(input_dict)
            )
        self.results['max_lambda'] = max_lambda_hists

        # calculate the TCPs
        tcp_methods = self.tcp_methods
        tcps = TransitionDictResults(
            {
                (trans.stateA, trans.stateB):
                tcp_methods[trans].from_ensemble_histograms(max_lambda_hists)
                for trans in tcp_methods
            },
            network=self.network
        )
        self.results['total_crossing_probability'] = tcps

        # calculate the CTPs
        ctps = self.ctp_method.from_weighted_trajectories(input_dict)
        self.results['conditional_transition_probability'] = ctps

        # calculate the transition probability from existing TCP, CTP
        #tp_m = self.transition_probability_methods
        #transition_probabilities = TransitionDictResults(
            #{
                #(t.stateA, t.stateB):
                #tp_m[t].from_intermediate_results(
                    #tcp=tcps[(t.stateA, t.stateB)],
                    #ctp=ctps
                #)
                #for t in tp_m
            #}
        #)
        #self.results['transition_probability'] = transition_probabilities

        pass

    def crossing_probability(self, ensemble):
        pass

    @property
    def conditional_transition_probability(self):
        pass

    @property
    def total_crossing_probability(self):
        pass
