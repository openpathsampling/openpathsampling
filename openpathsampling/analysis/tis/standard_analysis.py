import collections
import functools

import openpathsampling as paths
from openpathsampling.netcdfplus import StorableNamedObject
from openpathsampling.numerics import LookupFunction
import pandas as pd
import numpy as np

from .core import (MultiEnsembleSamplingAnalyzer, TransitionDictResults,
                   TISAnalysis, EnsembleHistogrammer)
from .crossing_probability import (
    FullHistogramMaxLambdas, TotalCrossingProbability
)
from .misc import ConditionalTransitionProbability
from .flux import MinusMoveFlux

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

class StandardTISAnalysis(TISAnalysis):
    """
    Standard TIS analysis: flux, TCP, CTP.

    This is what we call the "standard" TIS analysis. It splits the rate
    equation into a flux, a total crossing probability (calculated using
    ensemble crossing probability functions), and a conditional transition
    probability from the outermost interface.

    Whenever possible, this code allows you to use default values.

    For the **flux**, you must provide either a flux method or a move scheme
    (which will use the :class:`.MinusMoveFlux`).

    For the **conditional transition probability**, you may optionally
    provide a :class:`.ConditionalTransitionProbability` object, otherwise
    the code will create one for the outermost interfaces of each
    transition.

    For the **total crossing probability**, you must provide a dictionary
    for the ``max_lambda_calcs``. The keys of this dictionary are the
    sampling transitions; the values can either be an
    :class:`.EnsembleHistogrammer` (such as a
    :class:`.FullHistogramMaxLambdas`) or a dictionary of histogram
    parameters, in which case the histogram parameters will be passed to
    :class:`.FullHistogramMaxLambdas`. In addition, you may optionally
    provide a dictionary for ``combiners``, which maps the interface set
    within the sampling transitions to a combining function, such as
    :class:`.WHAM`. The default is to use :class:`.WHAM`.

    Parameters
    ----------
    network : :class:`.TISNetwork`
        the network to analyze
    steps : iterable of :class:`.MCStep`
        if given, the analysis is performed immediately using these steps;
        otherwise, the analysis can be performed later with
        :meth:`.calculate`
    flux_method : flux calculation method
        the method to use to calculate the flux; typical classes are
        :class:`.MinusMoveFlux` and :class:`.DictFlux`. Optional, but if not
        given then ``scheme`` must be given.
    scheme : :class:`.MoveScheme`
        used to create a :class:`.MinusMoveFlux` if ``flux_method`` is not
        provided. Not used if ``flux_method`` is given.
    ctp_method : :class:`ConditionalTransitionProbability`
        object for calculating the conditional transition probability
        (optional)
    max_lambda_calcs : dict
        determines how the ensemble crossing probability histograms are
        build. Keys are sampling transitions, and values can be either
        :class:`.EnsembleHistogrammer` subclasses, or a list of histogram
        parameters to pass to :class:`.FullHistogramMaxLambdas`.
    combiners : dict {:class:`.InterfaceSet`: combination method}
        links the interface set to the method that will be used to combine
        individual ensembles into the total crossing probability function.
        Default is to use :class:`.WHAM`.
    """
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

        self._progresser = paths.progress.SimpleProgress()

        if steps is not None:
            self.calculate(steps)

    def _set_progress(self, progress, leave):
        self._progresser.progress = progress  # set if string
        progress = self._progresser.progress  # get function
        prog_leave = functools.partial(progress, leave=True)
        prog_no_leave = functools.partial(progress, leave=False)
        if leave == 'all':
            self.flux_method.progress = prog_leave
            self.ctp_method.progress = prog_leave
            max_lambda_prog = prog_leave
        elif leave == 'default':
            self.flux_method.progress = prog_leave
            self.ctp_method.progress = prog_leave
            max_lambda_prog = prog_no_leave
        elif leave == 'none':
            self.flux_method.progress = prog_no_leave
            self.ctp_method.progress = prog_no_leave
            max_lambda_prog = prog_no_leave

        max_lambda_methods = [tcp.max_lambda_calc
                              for tcp in self.tcp_methods.values()]
        for max_lambda_m in max_lambda_methods:
            max_lambda_m.progress = max_lambda_prog

    @property
    def progress(self):
        if not hasattr(self._progresser, '_progress'):
            self.progress = 'default'
        return self._progresser.progress

    @progress.setter
    def progress(self, value):
        if value in ['all', 'default', 'none']:
            progress = 'tqdm'
            leave = value
        else:
            progress = value
            leave = 'default'
        self._set_progress(progress, leave)

    def from_weighted_trajectories(self, input_dict):
        """Calculate results from weighted trajectories dictionary.

        Parameters
        ----------
        input_dict : dict of {:class:`.Ensemble`: collections.Counter}
            ensemble as key, and a counter mapping each trajectory
            associated with that ensemble to its counter of time spent in
            the ensemble (output of `steps_to_weighted_trajectories`)

        Returns
        -------
        dict
            dictionary with all the results
        """
        # calculate the max_lambda hists
        max_lambda_calcs = [tcp_m.max_lambda_calc
                            for tcp_m in self.tcp_methods.values()]
        max_lambda_hists = {}
        label = "Crossing probability"
        for calc in self.progress(max_lambda_calcs, desc=label):
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
        """Crossing probability function for a given ensemble

        Parameters
        ----------
        ensemble : :class:`.Ensemble`
            the ensemble for which the crossing probability is desired

        Returns
        -------
        :class:`.LookupFunction`
            crossing probability function for the given ensemble
        """
        sampling_ens = self.network.sampling_ensemble_for[ensemble]
        all_max_lambdas = self._access_cached_result('max_lambda')
        max_lambda = all_max_lambdas[sampling_ens]
        return max_lambda.reverse_cumulative()


    @property
    def conditional_transition_probability(self):
        """``pandas.DataFrame``: conditional transition probabilities

        rows are ensemble names, columns are state names
        """
        ctp = self._access_cached_result('conditional_transition_probability')
        df = pd.DataFrame.from_dict(ctp, orient='index')
        df.index = [idx.name for idx in df.index]
        df.columns = [col.name for col in df.columns]
        return df

    @property
    def total_crossing_probability(self):
        """:class:`.LookupFunction`: total crossing probability
        """
        return self._access_cached_result('total_crossing_probability')
