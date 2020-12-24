import collections
import openpathsampling as paths
from openpathsampling.netcdfplus import StorableNamedObject
import pandas as pd
import numpy as np

from .core import MultiEnsembleSamplingAnalyzer

def flux_matrix_pd(flux_matrix, sort_method="default"):
    """Convert dict form of flux to a pandas.Series

    Parameters
    ----------
    flux_matrix : dict of {(state, interface): flux}
        the output of a flux calculation; flux out of state and through
        interface
    sort_method : callable or str
        method that takes a list of 2-tuple key from flux_matrix and returns
        a sorted list. Strings can be used to select internally-defined
        methods. Currently implemented: "default"
        (:meth:`.default_flux_sort`).

    Returns
    -------
    :class:`pandas.Series` :
        The flux represented in a pandas series
    """
    keys = list(flux_matrix.keys())
    known_method_names = {
        'default': default_flux_sort
    }
    if isinstance(sort_method, str):
        try:
            sort_method = known_method_names[sort_method.lower()]
        except KeyError:
            raise KeyError("Unknown sort_method name: " + str(sort_method))

    if sort_method is not None:
        ordered = sort_method(keys)
    else:
        ordered = keys
    values = [flux_matrix[k] for k in ordered]
    index_vals = [(k[0].name, k[1].name) for k in ordered]
    index = pd.MultiIndex.from_tuples(list(index_vals),
                                      names=["State", "Interface"])
    return pd.Series(values, index=index, name="Flux")


def default_flux_sort(tuple_list):
    """Default sort for flux pairs.

    Flux results are reported in terms of flux pairs like ``(state,
    interface)``. This sorts them using the ``.name`` strings for the
    volumes.
    """
    name_to_volumes = {(t[0].name, t[1].name): t for t in tuple_list}
    sorted_results = sorted(name_to_volumes.keys())
    return [name_to_volumes[key] for key in sorted_results]


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
            :meth:`.TrajectoryTransitionAnalysis.analyze_flux` (keys are
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
        flux_pairs = self.progress(self.flux_pairs, desc="Flux")
        for flux_pair in flux_pairs:
            (state, innermost) = flux_pair
            mover = flux_pair_to_minus_mover[flux_pair]
            calculator = transition_flux_calculators[flux_pair]
            minus_ens = flux_pair_to_minus_ensemble[flux_pair]
            # TODO: this won't work for SR minus, I don't think
            # (but neither would our old version)
            trajectories = [s.active[minus_ens].trajectory
                            for s in mover_to_steps[mover]]
            mover_trajs = self.progress(trajectories, leave=False)
            results[flux_pair] = calculator.analyze_flux(
                trajectories=mover_trajs,
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
            :meth:`.TrajectoryTransitionAnalysis.analyze_flux` (keys are
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
            :meth:`.TrajectoryTransitionAnalysis.analyze_flux` (keys are
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
            output of :meth:`.intermediates`

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
            output of :meth:`.intermediates`

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

