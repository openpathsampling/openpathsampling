import collections
import openpathsampling as paths
from openpathsampling.netcdfplus import StorableNamedObject

def steps_to_weighted_trajectories(steps, ensembles):
    """Bare function to convert to the weighted trajs dictionary.

    This prepares data for the faster analysis format. This preparation only
    need to be done once, and it will cover a lot of the analysis cases.
    """
    results = {e: collections.Counter() for e in ensembles}

    my_steps = steps
    # loop over blocks # TODO: add blocksize parameter, test various sizes
    block = collections.defaultdict(list)
    for step in my_steps:
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
        for ens in input_dict:
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
        if max_lambda_func is None:
            max_lambda_func = lambda t: max(transition.interfaces.cv(t))
            #max_lambda_func = transition.interfaces.max_cv  # TODO traj-cv
        self.lambdas = {e: l for (e, l) in zip(transition.ensembles,
                                               transition.lambdas)}
        super(FullHistogramMaxLambdas, self).__init__(
            ensembles=transition.ensembles,
            f=max_lambda_func,
            hist_parameters=hist_parameters
        )

#class PerEnsembleMaxLambdas(EnsembleHistogrammer):
    #def __init__(self, transition):
        #interfaces_lambdas = transition.interfaces.lambdas

class TotalCrossingProbability(MultiEnsembleSamplingAnalyzer):
    def __init__(self, max_lambda_calc, combiner=None):
        self.max_lambda_calc = max_lambda_calc
        if combiner is None:
            lambdas = self.max_lambda_calc.transition.interfaces.lambdas
            combiner = paths.numerics.WHAM(interfaces=lambdas)
        self.combiner = combiner

    def from_weighted_trajectories(self, input_dict):
        hists = self.max_lambda_calc.from_weighted_trajectories(input_dict)
        tcp_results = {}
        for trans in hists:
            df = paths.numerics.histograms_to_pandas_dataframe(
                hists.values(),
                fcn="reverse_cumulative"
            ).sort_index(axis=1)
            tcp = wham.wham_bam_histogram(df).to_dict()
            tcp_results[trans] = LookupFunction(tcp.keys(), tcp.values())
        return tcp_results


class ConditionalTransitionProbability(MultiEnsembleSamplingAnalyzer):
    def __init__(self, ensembles, states):
        self.states = states
        self.ensembles = ensembles

    def from_weighted_trajectories(self, input_dict):
        ctp = {}
        for ens in self.ensembles:
            acc = collecton.Counter()
            n_try = len(input_dict[ens])
            final_frames = [traj.get_as_proxy(1) for
                            traj in input_dict[ens].keys()]
            weights = input_dict[ens].values()
            for (f, w) in zip(final_frames, weights):
                acc += collection.Counter({s: w for s in states if s(f)})

            ctp[ens] = {s : float(acc[s]) / n_try}
            # TODO: add logging to report here
        return ctp


class TISAnalysis(StorableNamedObject):
    """
    Generic class for TIS analysis. One of these for each network.
    """
    def __init__(self, network, steps=None, flux=None, tcp=None, ctp=None):
        # NOTE: each of flux, ctp, tcp refer to the methods used; in
        # principle, these should have the option of being provided as a
        # single example (to be applied to all) or as a dict showing which
        # to apply to which analysis transition
        self.network = network

        # set default analysis behaviors
        if flux is None:
            flux_pairs = None  # TODO
            self.flux = MinusMoveFlux(flux_pairs)
        if tcp is None:
            self.tcp = {transition: PerEnsembleTCP(transition.ensembles)
                        for transition in self.network.sampling_transitions}
        if ctp is None:
            outermost_ensembles = [t.ensembles[-1]
                                   for t in self.network.sampling_transitions]
            self.ctp = ConditionalTransitionProbability(outermost_ensembles)

        self.transitions = network.transitions

        self.weighted_trajs = None  # default if there are no steps
        if steps is not None:
            self.weighted_trajs = steps_to_weighted_trajectories(
                steps,
                network.sampling_ensembles
            )
            pass  # do *all* the analysis!

    def rate_matrix(self, steps=None):
        pass

    def rate(self, from_state, to_state):
        pass

    def crossing_probability(self, ensemble):
        pass

    def flux(self, from_state, through_interface=None):
        pass

    def calc_flux(self, weighted_trajs):
        pass

    def calc_total_crossing_probability(self, weighted_trajs):
        sampling_tcps = {
            trans: self.tcp[traj].from_weighted_trajectories(weighted_trajs)
            for trans in self.network.sampling_transitions
        }

        for sampling_trans in sampling_tcps:
            pass

        # map the TCP from the sampling transitions to the analysis
        # transitions
        pass

    def calc_conditional_transition_probability(self, weighted_trajs):
        pass
