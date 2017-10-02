import collections
import openpathsampling as paths
from openpathsampling.netcdfplus import StorableNamedObject
from openpathsampling.numerics import LookupFunction
import pandas as pd
import numpy as np

from .core import MultiEnsembleSamplingAnalyzer, EnsembleHistogrammer


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
