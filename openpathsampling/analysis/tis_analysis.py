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

    Parameters
    ----------
    steps: iterable of :class:`.MCStep`
        steps to be analyzed
    ensembles: list of :class:`.Ensemble`
        ensembles to include in the list. Note: ensemble must be given!

    Returns
    -------
    dict of {:class:`.Ensemble`: collections.Counter}
        the result, with the ensemble as key, and a counter mapping each
        trajectory associated with that ensemble to its counter of time
        spent in the ensemble.
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
    """
    Abstract class for statistics from MC steps sampling multiple ensembles.

    Parameters
    ----------
    ensembles : list of :class:`.Ensemble`
        ensembles to be used in the calculation; can be overridden by
        :method:`.calculate`
    """
    def __init__(self, ensembles=None):
        self.ensembles = ensembles

    def calculate(self, steps, ensembles=None):
        """Perform the analysis, using `steps` as input.

        This is the main analysis for the abstract
        :class:`.MultiEnsembleSamplingAnalyzer`. Specific results depend on
        the specific subclass. Most objects simply need to override
        :method:`.from_weighted_trajectories` in order to obtain reasonable
        behavior.

        Parameters
        ----------
        steps : iterable of :class:`.MCStep`
            the steps to use as input for this analysis
        ensembles : list of :class:`.Ensemble
            ensembles to include in the calculation (other ensembles will be
            stripped); default is `None` meaning all ensembles given during
            initialization.

        Returns
        -------
        See .from_weighted_trajectories for this class.
        """
        if ensembles is None:
            ensembles = self.ensembles
        if ensembles is None:
            raise RuntimeError("If self.ensembles is not set, then "
                               + "ensembles must be given as argument to "
                               + "calculate")
        weighted_trajs = steps_to_weighted_trajectories(steps, ensembles)
        return self.from_weighted_trajectories(weighted_trajs)

    def from_weighted_trajectories(self, input_dict):
        """Calculate results from weighted trajectories dictionary.

        Must be implemented in subclass.
        """
        raise NotImplementedError

    @staticmethod
    def combine_results(result_1, result_2):
        """Combine two sets of results from this analysis.

        This can be used to combine results after parallelizing the
        analysis. The default is not implemented; it will only be
        implemented in cases where such a combination is feasible.
        """
        # to be used to simplify parallelization
        # TODO: implement this in subclasses in the future
        raise NotImplementedError

######## CALCULATING THE FLUX

class MinusMoveFlux(MultiEnsembleSamplingAnalyzer):
    """
    Calculating the flux from the minus move.

    Raises
    ------
    ValueError
        if the number of interface sets per minus move is greater than one.
        Cannot use Minus Move flux calculation with multiple interface set
        TIS.

    Parameters
    ----------
    scheme: :class:`.MoveScheme`
        move scheme that was used (includes information on the minus movers
        and on the network)
    flux_pairs: list of 2-tuple of :class:`.Volume`
        pairs of (state, interface) for calculating the flux out of the
        volume and through the state. Default is `None`, in which case the
        state and innermost interface are used.
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
        """
        Selects steps that used this object's minus movers
        """
        return [s for s in steps
                if s.change.canonical.mover in self.minus_movers
                and s.change.accepted]

    def trajectory_transition_flux_dict(self, minus_steps):
        """
        Main minus move-based flux analysis routine.

        Parameters
        ----------
        minus_steps: list of :class:`.MCStep`
            steps that used the minus movers

        Returns
        -------
        dict of {(:class:`.Volume, :class:`.Volume`): dict}
            keys are (state, interface); values are the result dict from
            :method:`.TrajectoryTransitionAnalysis.analyze_flux` (keys are
            strings 'in' and 'out', mapping to
            :class:`.TrajectorySegmentContainer` with appropriate frames.
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
        """Load from existing TrajectoryTransitionAnalysis calculations.

        Parameters
        ----------
        flux_dicts: dict of {(:class:`.Volume`, :class:`.Volume`): dict}
            keys are (state, interface); values are the result dict from
            :method:`.TrajectoryTransitionAnalysis.analyze_flux` (keys are
            strings 'in' and 'out', mapping to
            :class:`.TrajectorySegmentContainer` with appropriate frames.

        Returns
        -------
        dict of {(:class:`.Volume, :class:`.Volume`): float}
            keys are (state, interface); values are the associated flux
        """
        TTA = paths.TrajectoryTransitionAnalysis  # readability on 80 col
        return {k: TTA.flux_from_flux_dict(flux_dicts[k])
                for k in flux_dicts}

    def from_weighted_trajectories(self, input_dict):
        """Not implemented for flux calculation."""
        # this can't be done, e.g., in the case of the single replica minus
        # mover, where the minus trajectory isn't in the active samples
        raise NotImplementedError(
            "Can not calculate minus move from weighted trajectories."
        )

    def calculate(self, steps):
        """Perform the analysis, using `steps` as input.

        Parameters
        ----------
        steps : iterable of :class:`.MCStep`
            the steps to use as input for this analysis

        Returns
        -------
        dict of {(:class:`.Volume`, :class:`.Volume`): float}
            keys are (state, interface); values are the associated flux
        """
        intermediates = self.intermediates(steps)
        return self.calculate_from_intermediates(*intermediates)

    def intermediates(self, steps):
        """Calculate intermediates, using `steps` as input.

        Parameters
        ----------
        steps : iterable of :class:`.MCStep`
            the steps to use as input for this analysis

        Returns
        -------
        list (len 1) of dict of {(:class:`.Volume`, :class:`.Volume`): dict}
            keys are (state, interface); values are the result dict from
            :method:`.TrajectoryTransitionAnalysis.analyze_flux` (keys are
            strings 'in' and 'out', mapping to
            :class:`.TrajectorySegmentContainer` with appropriate frames.
        """
        minus_steps = self._get_minus_steps(steps)
        return [self.trajectory_transition_flux_dict(minus_steps)]

    def calculate_from_intermediates(self, *intermediates):
        """Perform the analysis, using intermediates as input.

        Parameters
        ----------
        intermediates :
            output of :method:`.intermediates`

        Returns
        -------
        dict of {(:class:`.Volume, :class:`.Volume`): float}
            keys are (state, interface); values are the associated flux
        """
        flux_dicts = intermediates[0]
        return self.from_trajectory_transition_flux_dict(flux_dicts)

class DictFlux(MultiEnsembleSamplingAnalyzer):
    """Pre-calculated flux, provided as a dict.

    Parameters
    ----------
    flux_dict: dict of {(:class:`.Volume`, :class:`.Volume`): float}
        keys are (state, interface) pairs; values are associated flux
    """
    def __init__(self, flux_dict):
        super(DictFlux, self).__init__()
        self.flux_dict = flux_dict

    def calculate(self, steps):
        """Perform the analysis, using `steps` as input.

        Parameters
        ----------
        steps : iterable of :class:`.MCStep`
            the steps to use as input for this analysis

        Returns
        -------
        dict of {(:class:`.Volume`, :class:`.Volume`): float}
            keys are (state, interface); values are the associated flux
        """
        return self.flux_dict

    def from_weighted_trajectories(self, input_dict):
        """Calculate results from weighted trajectories dictionary.

        For :class:`.DictFlux`, this ignores the input.

        Parameters
        ----------
        input_dict : dict of {:class:`.Ensemble`: collections.Counter}
            ensemble as key, and a counter mapping each trajectory
            associated with that ensemble to its counter of time spent in
            the ensemble.

        Returns
        -------
        dict of {(:class:`.Volume`, :class:`.Volume`): float}
            keys are (state, interface); values are the associated flux
        """
        return self.flux_dict

    def intermediates(self, steps):
        """Calculate intermediates, using `steps` as input.

        Parameters
        ----------
        steps : iterable of :class:`.MCStep`
            the steps to use as input for this analysis

        Returns
        -------
        list
            empty list; the method is a placeholder for this class
        """
        return []

    def calculate_from_intermediates(self, *intermediates):
        """Perform the analysis, using intermediates as input.

        Parameters
        ----------
        intermediates :
            output of :method:`.intermediates`

        Returns
        -------
        dict of {(:class:`.Volume, :class:`.Volume`): float}
            keys are (state, interface); values are the associated flux
        """
        return self.flux_dict

    @staticmethod
    def combine_results(result_1, result_2):
        """Combine two sets of results from this analysis.

        For :class:`.DictFlux`, the results must be identical.

        Parameters
        ----------
        result_1 : dict of {(:class:`.Volume, :class:`.Volume`): float}
            first set of results from a flux calculation
        result_2 : dict of {(:class:`.Volume, :class:`.Volume`): float}
            second set of results from a flux calculation

        Returns
        -------
        dict of {(:class:`.Volume, :class:`.Volume`): float}
            keys are (state, interface); values are the associated flux
        """
        if result_1 != result_2:
            raise RuntimeError("Combining results from different DictFlux")
        return result_1

########## GENERAL HISTOGRAMMING

class EnsembleHistogrammer(MultiEnsembleSamplingAnalyzer):
    """
    Generic code to calculate the properly weighted histograms of trajectory
    properties per ensemble.

    Parameters
    ----------
    ensembles: list of :class:`.Ensemble`
        ensembles to be included in the histogram
    f: callable
        the function to be histogrammed
    hist_parameters: dict
        allowed keys are 'bin_width' and 'bin_range'; value for 'bin_width'
        is a float; for 'bin_range' is a tuple with `(left_edge,
        right_edge)` (only left edge is used)
    """
    def __init__(self, ensembles, f, hist_parameters):
        super(EnsembleHistogrammer, self).__init__(ensembles)
        self.f = f
        self.hist_parameters = hist_parameters
        self.hists = {e: paths.numerics.Histogram(**self.hist_parameters)
                      for e in self.ensembles}

    def from_weighted_trajectories(self, input_dict):
        """Calculate results from a weighted trajectories dictionary.

        Parameters
        ----------
        input_dict : dict of {:class:`.Ensemble`: collections.Counter}
            ensemble as key, and a counter mapping each trajectory
            associated with that ensemble to its counter of time spent in
            the ensemble (output of `steps_to_weighted_trajectories`)

        Returns
        -------
        dict of {:class:`.Ensemble`: :class:`.numerics.Histogram`}
            calculated histogram for each ensemble
        """
        for ens in self.hists:
            trajs = list(input_dict[ens].keys())
            weights = list(input_dict[ens].values())
            data = [self.f(traj) for traj in trajs]
            self.hists[ens].histogram(data, weights)
        return self.hists


class PathLengthHistogrammer(EnsembleHistogrammer):
    """Histogramming path length distribution

    Parameters
    ----------
    ensembles: list of :class:`.Ensemble`
        ensembles to be included in the histogram
    hist_parameters: dict
        allowed keys are 'bin_width' and 'bin_range'; value for 'bin_width'
        is a float; for 'bin_range' is a tuple with `(left_edge,
        right_edge)` (only left edge is used)
    """
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

    This histograms the maximum value of lambda for each ensemble. One of
    these objects is made per transition.

    Parameters
    ----------
    transition: :class:`.TISTransition`
        the transition to be analyzed
    hist_parameters: dict
        Histogram parameters to use with this collective variable: allowed
        keys are 'bin_width' and 'bin_range'; value for 'bin_width' is a
        float; for 'bin_range' is a tuple with `(left_edge, right_edge)`
        (only left edge is used)
    max_lambda_func: callable
        function to use to map the trajectories to a histogram; default is
        `None`, which uses the maximum value of the order parameter
        associated with the interface set. Overriding this can be used if
        either (a) the interface set does not have an order parameter
        associated with it, or (b) you want to calculate the values along
        some other order parameter
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
    """
    Calculate the total crossing probability function.

    The total crossing probability function is generated by calculating the
    individual ensemble crossing probability functions (using, e.g.,
    :class:`.FullHistogramMaxLambdas`, and combining them using some
    combining method (default is :class:`.WHAM`). One of these objects is
    instantiated per transition.

    Parameters
    ----------
    max_lambda_calc: :class:`.EnsembleHistogrammer`
        usually :class:`.FullHistogramMaxLambdas`; object that creates the
        max lambda histograms for the ensembles associated with this
        transition.
    combiner: TODO
        class that combines multiple histograms (with restricted sampling)
        into a single result. If `None` (default), uses :class:`.WHAM`
    """
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
        """Calculate results from a weighted trajectories dictionary.

        Parameters
        ----------
        input_dict : dict of {:class:`.Ensemble`: collections.Counter}
            ensemble as key, and a counter mapping each trajectory
            associated with that ensemble to its counter of time spent in
            the ensemble (output of `steps_to_weighted_trajectories`)

        Results
        -------
        :class:`.LookupFunction`
            the total crossing probability function
        """
        hists = self.max_lambda_calc.from_weighted_trajectories(input_dict)
        return self.from_ensemble_histograms(hists)

    def from_ensemble_histograms(self, hists):
        """Calculate results from a dict of ensemble histograms.

        Parameters
        ----------
        hists : dict of {:class:`.Ensemble`: :class:`.numerics.Histogram`}
            histogram for each ensemble (from ``self.max_lambda_calc``)

        Returns
        -------
        :class:`.LookupFunction`
            the total crossing probability function
        """
        tcp_results = {}
        input_hists = [hists[ens] for ens in self.transition.ensembles]
        df = paths.numerics.histograms_to_pandas_dataframe(
            input_hists,
            fcn="reverse_cumulative"
        ).sort_index(axis=1)
        # TODO: remove WHAM-specific name here
        tcp = self.combiner.wham_bam_histogram(df).to_dict()
        return LookupFunction(tcp.keys(), tcp.values())


class ConditionalTransitionProbability(MultiEnsembleSamplingAnalyzer):
    """
    Calculate the conditional transition probability, P(B|A_m)

    The conditional transition probablity is the probability that paths from
    a given interface end in each possible state.

    Parameters
    ----------
    ensembles : list of :class:`.Ensemble`
        sampled ensembles to calculate the CTP for
    states : list of :class:`.Volume`
        possible final states for trajectories in the given ensembles
    """
    def __init__(self, ensembles, states):
        super(ConditionalTransitionProbability, self).__init__(ensembles)
        self.states = states

    def from_weighted_trajectories(self, input_dict):
        """Calculate results from a weighted trajectories dictionary.

        Parameters
        ----------
        input_dict : dict of {:class:`.Ensemble`: collections.Counter}
            ensemble as key, and a counter mapping each trajectory
            associated with that ensemble to its counter of time spent in
            the ensemble (output of `steps_to_weighted_trajectories`)

        Returns
        -------
        dict of {:class:`.Ensemble`: {:class:`.Volume`: float}}
            first key, an ensemble, selects the results from a given
            sampling ensemble; second key, a volume, selects the value for
            a given state. Value is the conditional transition probability
            for that state from that ensemble.
        """
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
    """
    Calculate the transition probability according to the TCP/CTP split.

    The transition probability is the probability that a path that starts
    from interface A_0 ends in some other state B, and is denoted P(B|A_0).
    What we call the "standard" approach to calculate this splits that
    probability into P(B|A_0) = P(B|A_m) P(A_m|A_0), where the first term in
    the product is the :class:`.ConditionalTransitionProbability` and the
    second term is determined by the :class:`.TotalCrossingProbability`. In
    practice, one instance of this object is created for each transition.

    Parameters
    ----------
    transition : :class:`.TISTransition`
        the transition to calculate the transition probability for
    tcp_method : :class:`.TotalCrossingProbability`
        the object to calculate the total crossing probability function;
        many details can be modified here including the details of the
        histograms that are generated or of the method used to combine the
        histograms
    ctp_method : :class:`.ConditionalTransitionProbability`
        the object to calculate the conditional transition probability
        (details on which ensembles and which states to analyze can be
        modified here)
    """
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
        """Calculate results from a weighted trajectories dictionary.

        Parameters
        ----------
        input_dict : dict of {:class:`.Ensemble`: collections.Counter}
            ensemble as key, and a counter mapping each trajectory
            associated with that ensemble to its counter of time spent in
            the ensemble (output of `steps_to_weighted_trajectories`)

        Returns
        -------
        float
            the transition probability for this transition
        """
        tcp = self.tcp_method.from_weighted_trajectories(input_dict)
        ctp = self.ctp_method.from_weighted_trajectories(input_dict)
        return self.from_intermediate_results(tcp, ctp)

    def from_intermediate_results(self, tcp, ctp):
        """Calculate results from intermediates.

        Parameters
        ----------
        tcp : :class:`.LookupFunction`
            results for the total crossing probability for this transition
        ctp : dict of {:class:`.Ensemble`: {:class:`.Volume`: float}}
            results for the conditional transition probability for this
            transition

        Returns
        -------
        float
            the transition probability for this transition
        """
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
    """Analysis result object for properties of a transition.

    Each value is associated with a specific (analysis/physical) transition.
    This object allows those values to be accessed either using (as a key,
    i.e., in square brackets) any of:

    * the transition object for the (initial_state, final_state) pair
    * the tuple of (initial_state, final_state) for the transition
    * the sampling transition object, if ``allow_sampling==True``; this is
      only desired if the quantity is only dependent on the sampling
      transition

    Note that results cannot be changed in this; a new object must be made
    if that is desired (but the original input can be accessed with the
    .results_dict attribute, and then modified as needed).

    Parameters
    ----------
    results_dict : dict of 2-tuple of :class:`.Volume` to float or Quantity
        dict connecting tuple of (initial_state, final_state) to the
        associated quantity
    network : :class:`.TransitionNetwork`
        the transition network associated with these results
    allow_sampling : bool
        whether to allow elements of network.sampling_transitions to be used
        as keys for to retrieve the stored results; this only makes sense if
        the stored result is only dependent on the sampling transition
    """
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
        """Output stored results as pandas.DataFrame

        Parameters
        ----------
        order : list of :class:`.Volume`
            order in which to list the states; if not used, the order may be
            unpredictable

        Returns
        -------
        :class:`pandas.DataFrame`
            DataFrame with initial states as rows and final states as
            columns
        """
        key_map = lambda key: key.name
        keys = list(self.results_dict.keys())
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

    In general, the TIS rate is split into the flux and the transition
    probability.

    Parameters
    ----------
    network : :class:`.TransitionNetwork`
        the reaction network to be analyzed
    flux_method : flux calculation method
        the method to use to calculate the flux; typical classes are
        :class:`.MinusMoveFlux` and :class:`.DictFlux`
    transition_probability_methods : dict of :class:`.Transition` to method
        the method for calculating the transition probability (one for each
        transition).
    """
    def __init__(self, network, flux_method, transition_probability_methods):
        self.network = network
        self.transitions = network.transitions
        self.flux_method = flux_method
        self.transition_probability_methods = transition_probability_methods
        self.results = {}

    def calculate(self, steps):
        """Perform the analysis, using `steps` as input.

        Parameters
        ----------
        steps : iterable of :class:`.MCStep`
            the steps to use as input for this analysis
        """
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
        """Calculate results from weighted trajectories dictionary.

        Parameters
        ----------
        input_dict : dict of {:class:`.Ensemble`: collections.Counter}
            ensemble as key, and a counter mapping each trajectory
            associated with that ensemble to its counter of time spent in
            the ensemble (output of `steps_to_weighted_trajectories`)
        """
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
        for (trans, transition_probability) in trans_prob.items():
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
        """dict of {(:class:`.Volume`, :class:`.Volume`): float}: keys are
        (state, interface); values are the associated flux
        """
        return self._access_cached_result('flux')

    def flux(self, from_state, through_interface=None):
        """Flux from a volume and through and interface.

        Shortcut to be used after the actual calculation has been performed.

        Parameters
        ----------
        from_state : :class:`.Volume`
            the volume the flux should start from
        through_interface : :class:`.Volume`
            the interface the flux should cross; default is None which uses
            the ``from_state`` volume

        Returns
        -------
        float or Quantity
            the flux out of the given state and through the given interface
        """
        fluxes = self._access_cached_result('flux')
        if through_interface is None:
            through_interface = from_state

        return fluxes[(from_state, through_interface)]

    def state_fluxes(self, from_state):
        """All fluxes associated with a given initial state.

        Shortcut to be used after the actual calculation has been performed.

        Parameters
        ----------
        from_state : :class:`.Volume`
            the volume the fluxes should start from

        Returns
        -------
        dict of 2-tuple of :class:`.Volume` to float
            dictionary of (state, interface) to the associated flux -- same
            as the flux dictionary given be :method:`.flux_matrix`, but only
            including the cases with the desired state volume
        """
        fluxes = self._access_cached_result('flux')
        state_keys = [k for k in fluxes.keys() if k[0] == from_state]
        return {k: fluxes[k] for k in state_keys}

    @property
    def transition_probability_matrix(self):
        """
        :class:`.TransitionDictResults`: matrix of transition probabilities
        """
        return self._access_cached_result('transition_probability')

    def transition_probability(self, from_state, to_state):
        """Transition probability between two states.

        Parameters
        ----------
        from_state : :class:`.Volume`
            initial state in the transition
        to_state : :class:`.Volume`
            final state in the transition

        Returns
        -------
        float
            transition probability for the `from_state`->`to_state`
            transition
        """
        trans_probs = self._access_cached_result('transition_probability')
        return trans_probs[(from_state, to_state)]

    def rate_matrix(self, steps=None):
        """Calculate the rate matrix.

        Parameters
        ----------
        steps : iterable of :class:`.MCStep`
            the steps from a simulation to use for calculating the rate. If
            `None` (default), then use the existing cached results.

        Returns
        -------
        :class:`.TransitionDictResults`
            the rate matrix
        """
        if steps is not None:
            self.calculate(steps)
        return self._access_cached_result('rate')

    def rate(self, from_state, to_state):
        """Rate for the transition between two states

        Parameters
        ----------
        from_state : :class:`.Volume`
            initial state in the transition
        to_state : :class:`.Volume`
            final state in the transition

        Returns
        -------
        float or Quantity
            rate for the `from_state`->`to_state` transition
        """
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
            raise RuntimeError("Must set max_lambda_calcs in "
                               + "StandardTISAnalysis")
        max_lambda_calc_dict = {}
        for (transition, calc) in max_lambda_calcs.items():
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
        for (trans, transition_probability) in trans_prob.items():
            trans_flux = fluxes[(trans.stateA, trans.interfaces[0])]
            rates[(trans.stateA, trans.stateB)] = \
                    trans_flux * transition_probability

        self.results['rate'] = TransitionDictResults(rates, self.network)
        return self.results


    def crossing_probability(self, ensemble):
        sampling_ens = self.network.sampling_ensemble_for[ensemble]
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
