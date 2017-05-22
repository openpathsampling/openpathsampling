import collections
import openpathsampling as paths
from openpathsampling.netcdfplus import StorableNamedObject
from openpathsampling.numerics import LookupFunction
import pandas as pd
import numpy as np

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
    def __init__(self, ensembles=None):
        self.ensembles = ensembles

    def calculate(self, steps, ensembles=None):
        if ensembles is None:
            ensembles = self.ensembles
        if ensembles is None:
            raise RuntimeError("If self.ensembles is not set, then "
                               + "ensembles must be given as argument to "
                               + "calculate")
        weighted_trajs = steps_to_weighted_trajectories(steps, ensembles)
        return self.from_weighted_trajectories(weighted_trajs)

    def from_weighted_trajectories(self, input_dict):
        raise NotImplementedError

    def combine_results(self, result_1, result_2):
        # to be used to simplify parallelization
        raise NotImplementedError

######## CALCULATING THE FLUX

class MinusMoveFlux(MultiEnsembleSamplingAnalyzer):
    """
    Parameters
    ----------
    network: :class:`.TransitionNetwork`
    flux_pairs: list of 2-tuple of :class:`.Volume`
    """
    def __init__(self, scheme, flux_pairs=None):
        super(MinusMoveFlux, self).__init__()
        # error string we'll re-use in a few places
        mistis_err_str = ("Cannot use minus move flux with multiple "
                          + "interface sets. ")
        self.scheme = scheme
        self.network = scheme.network
        self.minus_movers = scheme.movers['minus']
        for mover in self.minus_movers:
            n_innermost = len(mover.innermost_ensembles)
            if n_innermost != 1:
                raise ValueError(
                    mistis_err_str + "Mover " + str(mover) + " does not "
                    + "have exactly one innermost ensemble. Found "
                    + str(len(mover.innermost_ensembles)) + ")."
                )

        if flux_pairs is None:
            # get flux_pairs from network
            flux_pairs = []
            minus_ens_to_trans = self.network.special_ensembles['minus']
            for minus_ens in self.network.minus_ensembles:
                n_trans = len(minus_ens_to_trans[minus_ens])
                if n_trans > 1:  # pragma: no cover
                    # Should have been caught be the previous ValueError. If
                    # you hit this, something unexpected happened.
                    raise ValueError(mistis_err_str + "Ensemble "
                                     + repr(minus_ens) + " connects "
                                     + str(n_trans) + " transitions.")
                trans = minus_ens_to_trans[minus_ens][0]
                innermost = trans.interfaces[0]
                state = trans.stateA
                # a couple assertions as a sanity check
                assert minus_ens.state_vol == state
                assert minus_ens.innermost_vol == innermost
                flux_pairs.append((state, innermost))

        self.flux_pairs = flux_pairs

    def _get_minus_steps(self, steps):
        return [s for s in steps
                if s.change.canonical.mover in self.minus_movers
                and s.change.accepted]

    def trajectory_transition_flux_dict(self, minus_steps):
        """
        """
        # set up a few mappings that make it easier set up other things
        flux_pair_to_transition = {
            (trans.stateA, trans.interfaces[0]): trans
            for trans in self.network.sampling_transitions
        }

        flux_pair_to_minus_mover = {
            (m.minus_ensemble.state_vol, m.minus_ensemble.innermost_vol): m
            for m in self.minus_movers
        }

        minus_mover_to_flux_pair = {flux_pair_to_minus_mover[k]: k
                                    for k in flux_pair_to_minus_mover}

        flux_pair_to_minus_ensemble = {
            (minus_ens.state_vol, minus_ens.innermost_vol): minus_ens
            for minus_ens in self.network.minus_ensembles
        }

        # sanity checks -- only run once per analysis, so keep them in
        for pair in self.flux_pairs:
            assert pair in flux_pair_to_transition.keys()
            assert pair in flux_pair_to_minus_mover.keys()
        assert len(self.flux_pairs) == len(minus_mover_to_flux_pair)

        # organize the steps by mover used
        mover_to_steps = collections.defaultdict(list)
        for step in minus_steps:
            mover_to_steps[step.change.canonical.mover].append(step)

        # create the actual TrajectoryTransitionAnalysis objects to use
        transition_flux_calculators = {
            k: paths.TrajectoryTransitionAnalysis(
                transition=flux_pair_to_transition[k],
                dt=flux_pair_to_minus_mover[k].engine.snapshot_timestep
            )
            for k in self.flux_pairs
        }

        # do the analysis
        results = {}
        for flux_pair in self.flux_pairs:
            (state, innermost) = flux_pair
            mover = flux_pair_to_minus_mover[flux_pair]
            calculator = transition_flux_calculators[flux_pair]
            minus_ens = flux_pair_to_minus_ensemble[flux_pair]
            # TODO: this won't work for SR minus, I don't think
            # (but neither would our old version)
            trajectories = [s.active[minus_ens].trajectory
                            for s in mover_to_steps[mover]]
            results[flux_pair] = calculator.analyze_flux(
                trajectories=trajectories,
                state=state,
                interface=innermost
            )

        return results

    @staticmethod
    def from_trajectory_transition_flux_dict(flux_dicts):
        TTA = paths.TrajectoryTransitionAnalysis  # readability on 80 col
        return {k: TTA.flux_from_flux_dict(flux_dicts[k])
                for k in flux_dicts}

    def from_weighted_trajectories(self, input_dict):
        # this can't be done, e.g., in the case of the single replica minus
        # mover, where the final accepted trajectory
        raise NotImplementedError(
            "Can not calculate minus move from weighted trajectories."
        )

    def calculate(self, steps):
        intermediates = self.intermediates(steps)
        return self.calculate_from_intermediates(*intermediates)

    def intermediates(self, steps):
        minus_steps = self._get_minus_steps(steps)
        return [self.trajectory_transition_flux_dict(minus_steps)]

    def calculate_from_intermediates(self, *intermediates):
        flux_dicts = intermediates[0]
        return self.from_trajectory_transition_flux_dict(flux_dicts)

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

    def intermediates(self, steps):
        return []

    def calculate_from_intermediates(self, *intermediates):
        return self.flux_dict

    def combine_results(self, result_1, result_2):
        if result_1 != result_2:
            raise RuntimeError("Combining results from different DictFlux")
        return result_1

########## GENERAL HISTOGRAMMING

class EnsembleHistogrammer(MultiEnsembleSamplingAnalyzer):
    """
    Generic code to calculate the properly weighted histograms of trajectory
    properties per ensemble.
    """
    def __init__(self, ensembles, f, hist_parameters):
        super(EnsembleHistogrammer, self).__init__(ensembles)
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
            # TODO: check with PGB about these defaults
            hist_parameters = {'bin_width': 5, 'bin_range': (0, 1000)}

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
        transition = max_lambda_calc.transition
        super(TotalCrossingProbability, self).__init__(transition.ensembles)
        self.max_lambda_calc = max_lambda_calc
        self.transition = transition
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
        super(ConditionalTransitionProbability, self).__init__(ensembles)
        self.states = states

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
        self.ensembles = transition.ensembles
        self.outermost_ensemble = self.transition.ensembles[-1]
        interfaces = transition.interfaces
        self.outermost_lambda = interfaces.get_lambda(interfaces[-1])

    def from_weighted_trajectories(self, input_dict):
        tcp = self.tcp_method.from_weighted_trajectories(input_dict)
        ctp = self.ctp_method.from_weighted_trajectories(input_dict)
        return self.from_intermediate_results(tcp, ctp)

    def from_intermediate_results(self, tcp, ctp):
        outermost_ensemble_ctps = ctp[self.outermost_ensemble]
        try:
            outermost_ctp = outermost_ensemble_ctps[self.final_state]
        except KeyError:
            # no transition ends in that state
            outermost_ctp = 0.0  #  float('nan') # if you'd rather
        tcp_at_outermost = tcp(self.outermost_lambda)
        # TODO: log things here
        return outermost_ctp * tcp_at_outermost


class TransitionDictResults(StorableNamedObject):
    # allows you to use analysis transition, 2-tuple of states, or sampling
    # transition as the key to retrieve the stored results
    def __init__(self, results_dict, network, allow_sampling=True):
        # allow_sampling: can the sampling transitions be input?
        self.results_dict = results_dict
        self.network = network
        self.allow_sampling = allow_sampling

    def __iter__(self):
        return self.results_dict.__iter__()

    def __getitem__(self, key):
        if key in self.network.sampling_transitions and self.allow_sampling:
            key = self.network.sampling_to_analysis[key][0]
        try:
            key = (key.stateA, key.stateB)
        except AttributeError:
            # we have a stateA, stateB tuple
            pass
        return self.results_dict[key]

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

    def __str__(self):  # pragma: no cover
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
        flux_m = self.flux_method
        fluxes = flux_m.calculate(steps)
        self.results['flux'] = fluxes
        weighted_trajs = steps_to_weighted_trajectories(
            steps,
            self.network.sampling_ensembles
        )
        self.from_weighted_trajectories(weighted_trajs)

    def from_weighted_trajectories(self, input_dict):
        # dict of transition to transition probability
        tp_m = self.transition_probability_methods
        trans_prob = {t: tp_m[t].from_weighted_trajectories(input_dict)
                      for t in tp_m.keys()}
        self.results['transition_probability'] = TransitionDictResults(
            {(t.stateA, t.stateB) : trans_prob[t] for t in trans_prob},
            self.network
        )

        fluxes = self.flux_matrix
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
    def __init__(self, network, steps=None, flux_method=None, scheme=None,
                 ctp_method=None, max_lambda_calcs=None, combiners=None):
        # NOTE: each of flux, ctp, tcp refer to the methods used; in
        # principle, these should have the option of being provided as a
        # single example (to be applied to all) or as a dict showing which
        # to apply to which analysis transition

        # TODO: add logging to initialization to describe how the setup is
        # being interpreted

        # set default analysis behaviors
        if flux_method is None:
            if scheme is None:
                raise TypeError("StandardTISAnalysis requires either "
                                + "flux_method or scheme as argument.")
            flux_method = MinusMoveFlux(scheme)

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
            transition.interfaces: TotalCrossingProbability(
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
                tcp_method=self.tcp_methods[trans.interfaces],
                ctp_method=self.ctp_method
            )
            for trans in network.transitions.values()
        }

        super(StandardTISAnalysis, self).__init__(
            network=network,
            flux_method=flux_method,
            transition_probability_methods=trans_prob_methods
        )

        if steps is not None:
            self.calculate(steps)

    def from_weighted_trajectories(self, input_dict):
        # calculate the max_lambda hists
        max_lambda_calcs = [tcp_m.max_lambda_calc
                            for tcp_m in self.tcp_methods.values()]
        max_lambda_hists = {}
        for calc in max_lambda_calcs:
            calc_results = calc.from_weighted_trajectories(input_dict)
            # TODO: change this to a 2D mapping, CV and ensemble
            max_lambda_hists.update(calc_results)
        self.results['max_lambda'] = max_lambda_hists

        # calculate the TCPs
        # * raw_tcps take the sampling transitions, and map from interface set
        #   to results
        # * tcps create the map from (analysis) transition to (already
        #   calculated) results
        tcp_methods = self.tcp_methods
        raw_tcps = {
            ifaces:
            tcp_methods[ifaces].from_ensemble_histograms(max_lambda_hists)
            for ifaces in tcp_methods
        }
        tcps = TransitionDictResults(
            {(trans.stateA, trans.stateB): raw_tcps[trans.interfaces]
             for trans in self.network.transitions.values()},
            network=self.network
        )
        self.results['total_crossing_probability'] = tcps

        # calculate the CTPs
        ctps = self.ctp_method.from_weighted_trajectories(input_dict)
        self.results['conditional_transition_probability'] = ctps

        # calculate the transition probability from existing TCP, CTP
        fluxes = self.results['flux']
        tp_methods = self.transition_probability_methods
        trans_prob = {
            trans:
            tp_methods[trans].from_intermediate_results(
                tcp=tcps[(trans.stateA, trans.stateB)],
                ctp=ctps
            )
            for trans in tp_methods
        }
        transition_probabilities = TransitionDictResults(trans_prob,
                                                         self.network)
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


    def crossing_probability(self, ensemble, cv=None):
        sampling_ens = self.network.sampling_ensemble_for[ensemble]
        if cv is None:
            possible_cvs = [t.interfaces.cv
                            for t in self.network.sampling_transitions
                            if sampling_ens in t.ensembles]
            assert len(possible_cvs) == 1
            cv = possible_cvs[0]

        all_max_lambdas = self._access_cached_result('max_lambda')
        max_lambda = all_max_lambdas[sampling_ens]
        return max_lambda.reverse_cumulative()


    @property
    def conditional_transition_probability(self):
        ctp = self._access_cached_result('conditional_transition_probability')
        df = pd.DataFrame.from_dict(ctp, orient='index')
        df.index = [idx.name for idx in df.index]
        df.columns = [col.name for col in df.columns]
        return df

    @property
    def total_crossing_probability(self):
        return self._access_cached_result('total_crossing_probability')
