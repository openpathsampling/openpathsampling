import collections
import openpathsampling as paths
from openpathsampling.netcdfplus import StorableNamedObject

def steps_to_weighted_trajectories(steps, ensembles):
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
    def calculate(self, steps, ensembles):
        weighted_trajs = steps_to_weighted_trajectories(steps, ensembles)
        return self.from_weighted_trajectories(weighted_trajs)

    def from_weighted_trajectories(self, input_dict):
        raise NotImplementedError

#class MinusMoveFlux(MultiEnsembleSamplingAnalyzer):
    #def from_weighted_trajectories(self, input_dict):
        #pass

class DictFlux(MultiEnsembleSamplingAnalyzer):
    def __init__(self, flux_dict):
        super(DictFlux, self).__init__()
        self.flux_dict = flux_dict

    def calculate(self, steps):
        return self.flux_dict

    def from_weighted_trajectories(self, input_dict):
        return self.flux_dict

class EnsembleHistogrammer(MultiEnsembleSamplingAnalyzer):
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

class FullHistogramMaxLambas(EnsembleHistogrammer):
    def __init__(self, transition, hist_parameters, max_lambda_func=None):
        if max_lambda_func is None:
            max_lambda_func = lambda t: max(transition.interfaces.cv(t))
            #max_lambda_func = transition.interfaces.max_cv  # TODO traj-cv
        super(FullHistogramTCP, self).__init__(
            ensembles=transition.ensembles,
            f=max_lambda_func,
            hist_parameters=hist_parameters
        )


#class PerEnsembleMaxLambdas(EnsembleHistogrammer):
    #def __init__(self, transition):
        #interfaces_lambdas = transition.interfaces.lambdas

class PathLengthHistogrammer(EnsembleHistogrammer):
    def __init__(self, ensembles, hist_parameters=None):
        if hist_parameters is None:
            pass  # set defaults
        super(PathLengthHistogrammer, self).__init__(
            ensembles=ensembles,
            f=lambda t: len(t),
            hist_parameters=hist_parameters
        )

class ConditionalTransitionProbability(MultiEnsembleSamplingAnalyzer):
    def __init__(self, ensembles, states):
        pass

    def from_weighted_trajectories(self, input_dict):
        pass


class TISResults(StorableNamedObject):
    pass


class TISTransitionAnalysis(StorableNamedObject):
    def __init__(self, transition):
        pass


class TISAnalysis(StorableNamedObject):
    """
    Generic class for TIS analysis. One of these for each network.
    """
    def __init__(self, network, steps=None, flux=None, tcp=None, ctp=None):
        self.network = network
        if steps is not None:
            pass  # do *all* the analysis!

        # set default analysis behaviors
        if flux is None:
            flux_pairs = None  # TODO
            self.flux = MinusMoveFlux(flux_pairs)
        if tcp is None:
            self.tcp = {transition: PerEnsembleTCP(transition.ensembles)
                        for transition in self.network.sampling_transitions}
        if ctp is None:
            outermost_ensembles = [t.ensembles[-1]
                                   for t in self.sampling_transitions]
            self.ctp = ConditionalTransitionProbability(outermost_ensembles)

        self.transitions = network.transitions

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
        tcps = []
        for trans in self.network.sampling_transitions:
            tcp = self.tcp[trans].from_weighted_trajectories(weighted_trajs)
            tcps.append(tcp)

        pass

    def calc_conditional_transition_probability(self, weighted_trajs):
        pass
