import openpathsampling as paths
import numpy as np
from openpathsampling.high_level.move_strategy import levels
from openpathsampling.pathmover import SampleNaNError, SampleMaxLengthError
import openpathsampling.high_level.move_strategy as move_strategy
from functools import reduce


class SpringShootingSelector(paths.ShootingPointSelector):
    """
    Spring shooting point selector, symetric selector to be used in a
    spring shooting simulation. It uses a biased potential in the shape of
    min(1, e^(-k*i)) for a forward shooting move and min(1, e^(k*i)) for a
    backwards shooting move, where i is a frame number in the range
    [-delta_max, delta_max] and represents a shift (in frames) relative to
    the last accepted shooting frame index.

    Parameters
    ----------
    delta_max : int
        maximum shift in the shooting point index compared to the index of
        the last accepted shooting point
    k_spring : float
        k value to use in the biasing potential
    initial_guess : int
        Guess for the initial shooting point, if this is None it will default
        to len(trajectory)/2
    """
    def __init__(self, delta_max, k_spring, initial_guess=None):
        super(SpringShootingSelector, self).__init__()

        # Check if delta max is sane
        if delta_max <= 0:
            raise ValueError('delta_max should be bigger than 0')
        else:
            self.delta_max = delta_max

        # Check if k is sane
        if k_spring < 0:
            raise ValueError('k_spring should not be negative')
        else:
            self.k_spring = k_spring

        # Initiate the instance variable
        self.previous_snapshot = initial_guess
        self.trial_snapshot = initial_guess
        self.previous_trajectory = None

        # Make the bias lists
        self._fw_prob_list = self._spring_biases(self.delta_max,
                                                 -self.k_spring)
        self._bw_prob_list = self._spring_biases(self.delta_max,
                                                 self.k_spring)
        self._fw_total_bias = self.sum_spring_bias(self.delta_max,
                                                   -self.k_spring)
        self._bw_total_bias = self.sum_spring_bias(self.delta_max,
                                                   self.k_spring)

        # Check if the bias potentials are equal

        self._total_bias = self._fw_total_bias
        self.check_sanity()

    def f(self, snapshot, trajectory, direction=None):
        """
        Returns the unnormalized proposal probability of a snapshot

        Notes
        -----
        This needs a direction and only makes sense if snapshot is from the
        previous accepted trajectory.
        """
        if direction is None:
            raise NotImplementedError("f is not defined without a direction.")

        if str(direction).lower() not in {'forward', 'backward'}:
            raise NotImplementedError("direction must be either 'forward' or "
                                      "'backward'.")
        elif direction == 'forward':
            prob_list = self._fw_prob_list
        else:
            # Should be backward
            prob_list = self._bw_prob_list
        if trajectory is not self.previous_trajectory:
            raise NotImplementedError("f is not defined for any other "
                                      "trajectory than "
                                      "self.previous_trajectory.")
        if self.previous_snapshot is None:
            raise NotImplementedError("f is only defined if a previous index "
                                      "is known.")
        previous_shooting_index = self.previous_snapshot
        if previous_shooting_index < 0:
            previous_shooting_index += len(trajectory)

        idx = trajectory.index(snapshot)
        diff = (idx-previous_shooting_index)
        prob_idx = diff+self.delta_max
        if prob_idx < 0 or prob_idx >= len(prob_list):
            return 0
        else:
            return prob_list[prob_idx]

    def probability(self, snapshot, trajectory, direction=None):
        if direction is None:
            raise NotImplementedError("probability is not defined without a "
                                      "direction.")
        if str(direction).lower() not in {'forward', 'backward'}:
            raise NotImplementedError("direction must be either 'forward' or "
                                      "'backward'.")
        elif direction == 'forward':
            prob_total = self._fw_total_bias
        else:
            # Should be backwards
            prob_total = self._bw_total_bias
        return self.f(snapshot, trajectory, direction=direction)/prob_total

    def check_sanity(self):
        """
        Checks the sanity of the selector, making sure that de biases are
        symetric and the total sum is sane
        """

        if self._fw_prob_list != self._bw_prob_list[::-1]:
            raise RuntimeError("forward and backward biases are not equal")

        if sum(self._fw_prob_list[::-1]) != self._total_bias:
            raise RuntimeError("Sum of the biases changed")

    @staticmethod
    def _spring_biases(delta_max, k_spring):
        """
        Calculates the list of spring biases depending on delta_max and
        k_spring using the formula min(1, e^(k*i)) where i is in range
        [-delta_max,delta_max]

        Parameters
        ----------
        delta_max : int
            The maximum delta from 0
        k_spring : float
            The k value to use in the formula e^(k*i)

        Returns
        -------
        biases : list
            This is a list of biases of length 2*delta_max+1
        """
        return([min([1, np.exp(k_spring*float(i-delta_max))])
                for i in range(2*delta_max+1)])

    def sum_spring_bias(self, delta_max, k_spring):
        """
        Calculates the sum of the biases, given a delta_max and a k_spring.

        Parameters
        ----------
        delta_max : int
            The maximum delta from 0
        k_spring : int
            The k value to use in the formula e^(k*i)

        Returns
        -------
        sum_biases : float
            The sum of all the biases, summed from small to big
        """
        # Sum small to big to prevent float summing errors
        if k_spring < 0:
            return sum(self._spring_biases(delta_max, k_spring)[::-1])
        else:
            return sum(self._spring_biases(delta_max, k_spring))

    def probability_ratio(self, snapshot, initial_trajectory,
                          trial_trajectory):
        """
        Returns the acceptance probability of a trial trajectory, this is 1.0
        as long as a snapshot has been selected that is inside of the
        trajectory.

        Parameters
        ----------
        snapshot : :class:`.Snapshot`
            The shooting snapshot
        initial_trajectory : :class:`.Trajectory`
            The initial trajectory from which was shot
        trial_trajectory : :class:`.Trajectory`
            The generated trial trajectory

        Returns
        -------
            float:
                1.0 if anacceptable snapshot has been chosen 0.0 otherwise
        """
        # Check if an acceptable snapshot was selected
        if self._acceptable_snapshot:
            return 1.0
        else:
            return 0.0

    def pick(self, trajectory, direction=None):
        """
        Picks a frame index to simulate from

        Parameters
        ----------
        trajectory : :class:`.Trajectory`
            The trajectory from which a frame has to be picked
        direction : string
            The direction of the engine mover, either forward or backward

        Returns
        -------
        index : int
            The index of the selected shooting snapshot in the given trajectory
        """
        self.check_sanity()
        # Get delta_max
        delta_max = self.delta_max

        # Check if last trajectory was accepted and update variable
        if trajectory != self.previous_trajectory:
            self.previous_snapshot = self.trial_snapshot

        # Check if previous snapshot is present, default to int devision
        # otherwise
        if self.previous_snapshot is None:
            previous_shooting_index = int(len(trajectory)/2)
        else:
            previous_shooting_index = self.previous_snapshot

        # A negative shooting index is needed to keep track of the right index
        # if the backward part of a trajectory changed length
        if previous_shooting_index < 0:
            previous_shooting_index += len(trajectory)

        # Select the right probability list to use
        if direction == 'forward':
            prob_list = self._fw_prob_list
        elif direction == 'backward':
            prob_list = self._bw_prob_list
        else:
            raise RuntimeError('Pick called without direction')

        # Select the next index from the probability list and shift to range
        # from -delta_max to delta_max, 0 being the shooting index of the last
        # accepted trajectory
        rand = self._rng.random() * self._total_bias
        dframe = 0
        prob = prob_list[0]
        while prob <= rand and dframe < len(prob_list):
            dframe += 1
            prob += prob_list[dframe]

        index = previous_shooting_index + (dframe - delta_max)

        # Check if a valid index has been selected (exclude state frames)
        if index <= 0 or index >= len(trajectory)-1:
            # This will always result in a 0 md-step try
            # either start at 0 for backwards or the last frame for forward
            if direction == "backward":
                index = 0
            else:
                index = len(trajectory)-1
            self._acceptable_snapshot = False
            # Needed to prevent selecting an out of range index
        else:
            self._acceptable_snapshot = True

        # Keeps track of the right index if the backward part changes length
        if direction == "backward":
            self.trial_snapshot = index - len(trajectory)
        else:
            self.trial_snapshot = index

        # Keep track of the current trajectory
        self.previous_trajectory = trajectory
        return index

    def restart_from_step(self, step):
        """
        This restarts the selector, based on the move details of the given
        step.

        Properties
        ----------
        step : :class:`.MCStep`
            The Monte Carlo step to restart the simulation from
        """
        # Get the details
        details = step.change.canonical.details

        # Try to reset instance variables
        try:
            self.previous_trajectory = details.initial_trajectory
            self.previous_snapshot = details.last_accepted_shooting_index
            self.trial_snapshot = details.shooting_index
            if details.direction == 'backward':
                self.trial_snapshot -= len(self.previous_trajectory)
        except AttributeError:
            raise RuntimeError("Tried to restart from a step that was not made\
                                by the spring shooting algorithm")


class SpringMover(paths.pathmover.EngineMover):
    """A shooting sample generator for the spring shooting algorithm

    Properties
    ----------
    ensemble : :class:`.Ensemble`
        Ensemble to simulate in
    selector : :class:`.ShootingPointSelector`
        The selector to select the shooting points
    engine : :class:`.DynamicsEngine`
        The dynamics engine to run the dynamics
    """

    default_engine = None
    reject_max_length = True

    def __init__(self, ensemble, selector=None, engine=None):
        super(SpringMover, self).__init__(
            ensemble=ensemble,
            target_ensemble=ensemble,
            selector=selector,
            engine=engine
        )

    @property
    def direction(self):
        raise NotImplementedError('SpringMover should have a direction')

    def __call__(self, input_sample):

        # Get the input trajectory and pick a frame
        initial_trajectory = input_sample.trajectory
        shooting_index = self.selector.pick(initial_trajectory,
                                            direction=self.direction)
        # Convert internal counter to a positive number for details
        prev_index = self.selector.previous_snapshot
        if prev_index is not None and prev_index < 0:
            prev_index += len(initial_trajectory)

        # Reject imposible snapshots before running any dynamics
        if not self.selector._acceptable_snapshot:
            trial_trajectory = initial_trajectory
            trial, details = self._build_sample(input_sample,
                                                shooting_index,
                                                trial_trajectory,
                                                'invalid_index')
            trials = [trial]
            details.update({'rejection_reason': 'invalid_index',
                            'shooting_index': shooting_index,
                            'last_accepted_shooting_index': prev_index,
                            'direction': self.direction
                            })
            return trials, details

        # Run dynamics and catch the engine errors
        try:
            trial_trajectory, run_details = self._run(initial_trajectory,
                                                      shooting_index)

        except paths.engines.EngineNaNError as err:
            trial, details = self._build_sample(
                input_sample, shooting_index, err.last_trajectory, 'nan')

            raise SampleNaNError('Sample with NaN', trial, details)

        except paths.engines.EngineMaxLengthError as err:
            trial, details = self._build_sample(
                input_sample, shooting_index, err.last_trajectory,
                'max_length')

            if self.reject_max_length:
                raise SampleMaxLengthError('Sample with MaxLength',
                                           trial, details)

        # Make the sample
        else:
            trial, details = self._build_sample(
                input_sample, shooting_index, trial_trajectory)

        # Make a list of trials and update the details.
        trials = [trial]
        details.update(run_details)
        details.update({'shooting_index': shooting_index,
                        'last_accepted_shooting_index': prev_index,
                        'direction': self.direction
                        })

        return trials, details


class ForwardSpringMover(SpringMover):
    """A forward shooting sample generator for the spring shooting algorithm

    Properties
    ----------
    ensemble : :class:`.Ensemble`
        Ensemble to simulate in
    selector : :class:`.ShootingPointSelector`
        The selector to select the shooting points
    engine : :class:`.DynamicsEngine`
        The dynamics engine to run the dynamics
    """

    default_engine = None
    reject_max_length = True

    def __init__(self, ensemble, engine=None, selector=None):
        super(ForwardSpringMover, self).__init__(
            ensemble=ensemble,
            selector=selector,
            engine=engine
        )

    @property
    def direction(self):
        return 'forward'


class BackwardSpringMover(SpringMover):
    """A backward shooting sample generator for the spring shooting algorithm

    Properties
    ----------
    ensemble : :class:`.Ensemble`
        Ensemble to simulate in
    selector : :class:`.ShootingPointSelector`
        The selector to select the shooting points
    engine : :class:`.DynamicsEngine`
        The dynamics engine to run the dynamics
    """

    default_engine = None
    reject_max_length = True

    def __init__(self, ensemble, engine=None, selector=None):
        super(BackwardSpringMover, self).__init__(
            ensemble=ensemble,
            selector=selector,
            engine=engine
        )

    @property
    def direction(self):
        return 'backward'


class SpringShootingMover(paths.pathmover.SpecializedRandomChoiceMover):
    """
    A Sample mover implementing the spring shooting algorithm

    Parameters
    ----------
    ensemble : :class:`.Ensemble`
        The simulation ensemble
    delta_max : int
        maximum shift in the shooting point index compared to the index of
        the last accepted shooting point
    k_spring : float
        k value to use in the biasing potential
    engine : :class:`.DynamicsEngine`
        The dynamics engine to run the dynamics
    initial_guess : int
        Guess for the initial shooting point, if this is None it will default
        to len(trajectory)/2
    """

    def __init__(self, ensemble, delta_max, k_spring, engine=None,
                 initial_guess=None):
        # Initialize the selector beforehand, so that both movers use the same
        # instance of the selector and it can correctly keep track of the
        # previous selected shots, independent of the mover used
        selector = SpringShootingSelector(delta_max=delta_max,
                                          k_spring=k_spring,
                                          initial_guess=initial_guess)

        movers = [
            ForwardSpringMover(
                ensemble=ensemble,
                engine=engine,
                selector=selector
            ),
            BackwardSpringMover(
                ensemble=ensemble,
                engine=engine,
                selector=selector
            )
        ]
        super(SpringShootingMover, self).__init__(
            movers=movers
        )

    @property
    def ensemble(self):
        return self.movers[0].ensemble

    @property
    def selector(self):
        return self.movers[0].selector


class SpringShootingStrategy(move_strategy.SingleEnsembleMoveStrategy):

    """
    Strategy for SpringShooting. Using the spring shooting point selector.
    Parameters
    ----------
    ensembles : list of :class:`.Ensemble`
        ensembles for which this strategy applies; None gives default
        behavior
    engine : :class:`.DynamicsEngine`
        engine for the dynamics
    group : str
        mover group name, default "shooting"
    replace : bool
        whether to replace existing movers in the group; default True
    """
    _level = levels.MOVER

    def __init__(self, delta_max, k_spring, initial_guess=None,
                 ensembles=None, engine=None, group="shooting", replace=True):
        super(SpringShootingStrategy, self).__init__(
            ensembles=ensembles, group=group, replace=replace
        )
        self.engine = engine
        self.delta_max = delta_max
        self.k_spring = k_spring
        self.initial_guess = initial_guess

    def make_movers(self, scheme):
        ensemble_list = self.get_init_ensembles(scheme)
        ensembles = reduce(list.__add__, map(lambda x: list(x), ensemble_list))
        movers = []
        for ens in ensembles:
            mover = SpringShootingMover(
                ensemble=ens,
                delta_max=self.delta_max,
                k_spring=self.k_spring,
                initial_guess=self.initial_guess,
                engine=self.engine
            )
            mover.named("SpringShootingMover " + str(ens.name))
            movers.append(mover)

        return movers


class SpringShootingMoveScheme(paths.high_level.move_scheme.MoveScheme):
    """
    A move scheme implementing the spring shooting algorithm for one way TPS

    Parameters
    ----------
    network : :class:`.TransitionNetworks`
    ensemble : :class:`.Ensemble`
        The simulation ensemble
    delta_max : int
        maximum shift in the shooting point index compared to the index of
        the last accepted shooting point
    k_spring : float
        k value to use in the biasing potential
    initial_guess : int
        Guess for the initial shooting point, if this is None it will default
        to len(trajectory)/2
    ensembles : list of :class:`.Ensemble`
        ensembles for which this move scheme applies. None gives default
        behaviour.
    engine : :class:`.DynamicsEngine`
        The dynamics engine to run the dynamics
    """
    def __init__(self, network, delta_max, k_spring, initial_guess=None,
                 ensembles=None, engine=None):
        self.delta_max = delta_max
        self.k_spring = k_spring
        self.initial_guess = initial_guess
        super(SpringShootingMoveScheme, self).__init__(network)
        self.append(SpringShootingStrategy(ensembles=ensembles,
                                           delta_max=delta_max,
                                           k_spring=k_spring,
                                           initial_guess=initial_guess,
                                           engine=engine))
        self.append(move_strategy.OrganizeByMoveGroupStrategy())

    def to_dict(self):
        ret_dict = {
            'movers': self.movers,
            'network': self.network,
            'choice_probability': self.choice_probability,
            'real_choice_probability': self.real_choice_probability,
            'balance_partners': self.balance_partners,
            'root_mover': self.root_mover,
            'delta_max': self.delta_max,
            'k_spring': self.k_spring,
            'initial_guess': self.initial_guess
        }
        return ret_dict

    @classmethod
    def from_dict(cls, dct):
        scheme = cls.__new__(cls)
        # noinspection PyArgumentList
        scheme.__init__(network=dct['network'],
                        delta_max=dct['delta_max'],
                        k_spring=dct['k_spring'],
                        initial_guess=dct['initial_guess'])
        scheme.movers = dct['movers']
        scheme.choice_probability = dct['choice_probability']
        scheme._real_choice_probability = dct['real_choice_probability']
        scheme.balance_partners = dct['balance_partners']
        scheme.root_mover = dct['root_mover']
        return scheme
