import collections
import openpathsampling as paths
from openpathsampling.netcdfplus import StorableNamedObject
from openpathsampling.progress import SimpleProgress
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
            result.at[key_map(k[0]), key_map(k[1])] = self.results_dict[k]

        return result

    def __str__(self):  # pragma: no cover
        return self.to_pandas().__str__()

    def __repr__(self):
        return self.to_pandas().__repr__()

class MultiEnsembleSamplingAnalyzer(SimpleProgress, StorableNamedObject):
    """
    Abstract class for statistics from MC steps sampling multiple ensembles.

    Parameters
    ----------
    ensembles : list of :class:`.Ensemble`
        ensembles to be used in the calculation; can be overridden by
        :meth:`.calculate`
    """
    def __init__(self, ensembles=None):
        self.ensembles = ensembles

    def calculate(self, steps, ensembles=None):
        """Perform the analysis, using `steps` as input.

        This is the main analysis for the abstract
        :class:`.MultiEnsembleSamplingAnalyzer`. Specific results depend on
        the specific subclass. Most objects simply need to override
        :meth:`.from_weighted_trajectories` in order to obtain reasonable
        behavior.

        Parameters
        ----------
        steps : iterable of :class:`.MCStep`
            the steps to use as input for this analysis
        ensembles : list of :class:`.Ensemble`
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
        steps = self.progress(steps, desc="Weighted trajectories")
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
    _label = "Ensembles"
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
        hists = self.progress(self.hists, desc=self._label)
        for ens in hists:
            trajs = input_dict[ens].keys()
            weights = list(input_dict[ens].values())
            data = [self.f(traj)
                    for traj in self.progress(trajs, leave=False)]
            self.hists[ens].histogram(data, weights)
        return self.hists


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
            as the flux dictionary given be :meth:`.flux_matrix`, but only
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
