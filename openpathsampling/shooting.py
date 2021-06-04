import math
import logging

import numpy as np

from openpathsampling.netcdfplus import StorableNamedObject
from openpathsampling import default_rng

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


class ShootingPointSelector(StorableNamedObject):
    def __init__(self):
        # Assign rng, so it can be set to something else
        self._rng = default_rng()
        super(ShootingPointSelector, self).__init__()

    def f(self, snapshot, trajectory):
        """
        Returns the unnormalized proposal probability of a snapshot

        Notes
        -----
        In principle this is an collectivevariable so we could easily add
        caching if useful
        """
        return 1.0

    def probability(self, snapshot, trajectory):
        sum_bias = self.sum_bias(trajectory)
        if sum_bias > 0.0:
            return self.f(snapshot, trajectory) / sum_bias
        else:
            return 0.0

    def probability_ratio(self, snapshot, old_trajectory, new_trajectory):
        p_old = self.probability(snapshot, old_trajectory)
        p_new = self.probability(snapshot, new_trajectory)
        return p_new / p_old

    def _biases(self, trajectory):
        """
        Returns a list of unnormalized proposal probabilities for all
        snapshots in trajectory
        """
        return [self.f(s, trajectory) for s in trajectory]

    def sum_bias(self, trajectory):
        """
        Returns the unnormalized probability probability of a trajectory.
        This is just the sum of all proposal probabilities in a trajectory.

        Notes
        -----
        For a uniform distribution this is proportional to the length of the
        trajectory. In this case we can estimate the maximal accepted
        trajectory length for a given acceptance probability.

        After we have generated a new trajectory the acceptance probability
        only for the non-symmetric proposal of different snapshots is given
        by `probability(old_trajectory) / probability(new_trajectory)`
        """

        return sum(self._biases(trajectory))

    def pick(self, trajectory):
        """
        Returns the index of the chosen snapshot within `trajectory`

        Notes
        -----
        The native implementation is very slow. Simple picking algorithm
        should override this function.
        """

        prob_list = self._biases(trajectory)
        sum_bias = sum(prob_list)

        rand = self._rng.random() * sum_bias
        idx = 0
        prob = prob_list[0]
        while prob <= rand and idx < len(prob_list):
            idx += 1
            prob += prob_list[idx]

        return idx


class GaussianBiasSelector(ShootingPointSelector):
    r"""
    A selector that biases according to a Gaussian along specified
    :class:`.CollectiveVariable`, with mean ``l_0`` and width parameter
    ``alpha``.  That is, for snapshot :math:`x` and CV :math:`\lambda`, the
    selection probability for each frame is weighted according to the
    function

    .. math::

        P_\text{sel}(x) \propto \exp(-\alpha (\lambda(x) - l_0)^2)

    Note that normalization here depends on the trajectory that the
    snapshot is a part of: the sum of the probabilities for all frames
    is 1, which gives a different normalization constant than the standard
    Gaussian distribution normalization, and exact probabilities for
    selecting a given snapshot will change depending on the trajectory it is
    a part of.

    Parameters
    ----------
    collectivevariable : :class:`.CollectiveVariable`
        the axis to use for the Gaussian
    alpha : float
        the width of the Gaussian
    l_0 : float
        the center of the Gaussian
    """
    def __init__(self, collectivevariable, alpha=1.0, l_0=0.5):
        super(GaussianBiasSelector, self).__init__()
        self.collectivevariable = collectivevariable
        self.alpha = alpha
        self.l_0 = l_0

    def f(self, snapshot, trajectory):
        l_s = self.collectivevariable(snapshot)
        return math.exp(-self.alpha * (l_s - self.l_0) ** 2)


class BiasedSelector(ShootingPointSelector):
    """General biased shooting point selector

    Takes any function (wrapped in an OPS CV) and uses that as the bias for
    selecting the shooting point.

    Parameters
    ----------
    func : :class:`.CollectiveVariable`
        A function wrapped in an OPS CV which gives the relative bias.
    """
    def __init__(self, func):
        super(BiasedSelector, self).__init__()
        self.func = func

    def f(self, snapshot, trajectory):
        return self.func(snapshot)


class UniformSelector(ShootingPointSelector):
    """
    Selects random frame in range `pad_start` to `len(trajectory-pad_end`.

    Attributes
    ----------
    pad_start : int
        number of frames at beginning of trajectory to be excluded from
        selection
    pad_end : int
        number of frames at end of trajectory to be excluded from selection
    """

    def __init__(self, pad_start=1, pad_end=1):
        super(UniformSelector, self).__init__()
        self.pad_start = pad_start
        self.pad_end = pad_end

    def f(self, frame, trajectory=None):
        return 1.0

    def sum_bias(self, trajectory):
        return float(len(trajectory) - self.pad_start - self.pad_end)

    def pick(self, trajectory):
        idx = self._rng.integers(self.pad_start,
                                len(trajectory) - self.pad_end)
        return idx


class InterfaceConstrainedSelector(ShootingPointSelector):
    """
    Selects first frame outside of volume.

    Parameters
    ----------
    volume : :class:`.Volume`
        defines Volume for which the first frame outside of this interface
        volume is found
    """

    def __init__(self, volume):
        super(InterfaceConstrainedSelector, self).__init__()
        self.volume = volume

    def f(self, frame, trajectory=None):
        idx = trajectory.index(frame)
        if idx == self.pick(trajectory):
            return 1.0
        else:
            return 0.0

    def sum_bias(self, trajectory):
        return 1.0

    def pick(self, trajectory):
        for idx, frame in enumerate(trajectory):
            if not self.volume(frame):
                break
        if idx == len(trajectory)-1 and self.volume(frame):
            raise RuntimeError("Interface constrained shooting move did "
                               " not find valid crossing point")

        return idx


class FinalFrameSelector(ShootingPointSelector):
    """
    Pick final trajectory frame as shooting point.

    This is used for "forward" extension in, e.g., the minus move.
    """
    def f(self, frame, trajectory):
        if trajectory.index(frame) == len(trajectory) - 1:
            return 1.0
        else:
            return 0.0

    def pick(self, trajectory):
        return len(trajectory)-1

    def probability(self, snapshot, trajectory):  # pragma: no cover
        return 1.0  # there's only one choice

    def probability_ratio(self, snapshot, old_trajectory, new_trajectory):
        # must be matched by a final-frame selector somewhere
        return 1.0


class FirstFrameSelector(ShootingPointSelector):
    """
    Pick first trajectory frame as shooting point.

    This is used for "backward" extension in, e.g., the minus move.
    """

    def f(self, frame, trajectory):
        if trajectory.index(frame) == 0:
            return 1.0
        else:
            return 0.0

    def pick(self, trajectory):
        return 0

    def probability(self, snapshot, trajectory):  # pragma: no cover
        return 1.0  # there's only one choice

    def probability_ratio(self, snapshot, old_trajectory, new_trajectory):
        # must be matched by a first-frame selector somewhere
        return 1.0
