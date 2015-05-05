"""Class to have trial generating instances

"""
__author__ = 'jan-hendrikprinz'

import openpathsampling as paths
import logging
import random
from openpathsampling.todict import ops_object

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')

class SampleGenerator(paths.PathMover):
    engine = None

    @staticmethod
    def MCAcceptor(trials, ensembles):
        trial_dict = dict()
        for trial in trials:
            trial_dict[trial.ensemble] = trial

        accepted = True
        probability = 1.0

        for ens in ensembles:
            if ens in trial_dict:
                sample = trial_dict[ens]
            else:
                accepted = False
                break

            valid = ens(sample)
            if not valid:
                # one sample not valid reject
                accepted = False
                break
            else:
                probability *= sample.bias

        rand = random.random()

        if rand > probability:
            # rejected
            accepted = False

        details = paths.MoveDetails(
            total_acceptance = probability,
            random_value = rand
        )

        return accepted, details

    def __init__(self, in_ensembles, out_ensembles):
        super(SampleGenerator, self).__init__(in_ensembles)
        self.out_ensembles = out_ensembles
        self.acceptor = SampleGenerator.MCAcceptor

    def move(self, globalstate):
        # we use self.ensembles to pick the samples we want to use
        samples = [ globalstate[ens] for ens in self.ensembles ]

        # pass these samples to the generator
        trials = self._generate(*samples)

        # accept/reject
        accepted, details = self._accept(trials)

        # and return a PMC
        if accepted:
            return paths.AcceptedSamplePathMoveChange(
                mover = self,
                details = details
            )
        else:
            return paths.RejectedSamplePathMoveChange(
                mover = self,
                details = details
            )

    def __call__(self, *args):
        return args

    def _accept(self, trials):
        return self.acceptor(trials, self.out_ensembles)

###############################################################################
# SHOOTERS
###############################################################################

class ShootingGenerator(SampleGenerator):
    def __init__(self, selector, ensemble):
        super(ShootingGenerator, self).__init__(
            in_ensembles=[ensemble],
            out_ensembles=[ensemble]
        )
        self.selector = selector

    def __call__(self, trial):
        initial_trajectory = trial.trajectory

        dynamics_ensemble = trial.ensemble
        replica = trial.replica

        initial_point = self.selector.pick(initial_trajectory)
        trial_point = self._shoot(initial_point)

        bias = initial_point.sum_bias / trial_point.sum_bias

        trial_details = paths.SampleDetails(
            initial_point=initial_point,
            trial_point=trial_point,
        )

        trial = paths.Sample(
            replica=replica,
            trajectory=trial_point.trajectory,
            ensemble=dynamics_ensemble,
            parent=trial,
            details=trial_details,
            mover=self,
            bias=bias
        )

        trials = [trial]

        return trials

    def _shoot(self, shooting_point, ensemble):
        return shooting_point


class ForwardShootGenerator(ShootingGenerator):
    def _shoot(self, shooting_point, ensemble):
        shoot_str = "Shooting {sh_dir} from frame {fnum} in [0:{maxt}]"
        logger.info(shoot_str.format(fnum=shooting_point.index,
                                     maxt=len(shooting_point.trajectory)-1,
                                     sh_dir="forward",
                                    ))

        # Run until one of the stoppers is triggered
        partial_trajectory = SampleGenerator.engine.generate(
            shooting_point.initial_point.copy(),
            running = [
                paths.ForwardAppendedTrajectoryEnsemble(
                    ensemble,
                    shooting_point.trajectory[0:shooting_point.index]
                ).can_append,
                self._length_stopper.can_append
            ]
        )

        trial_trajectory = shooting_point.trajectory[0:shooting_point.index] + partial_trajectory

        trial_point = paths.ShootingPoint(
            shooting_point.selector,
            trial_trajectory,
            shooting_point.index
        )

        return trial_point


class BackwardShootingGenerator(ShootingGenerator):
    def _shoot(self, shooting_point, ensemble):
        shoot_str = "Shooting {sh_dir} from frame {fnum} in [0:{maxt}]"
        logger.info(shoot_str.format(fnum=shooting_point.index,
                                     maxt=len(shooting_point.trajectory)-1,
                                     sh_dir="backward",
                                    ))

        # Run until one of the stoppers is triggered
        partial_trajectory = SampleGenerator.engine.generate(
            shooting_point.initial_point.reversed_copy(),
            running = [
                paths.BackwardPrependedTrajectoryEnsemble(
                    ensemble,
                    shooting_point.trajectory[shooting_point.index + 1:]
                ).can_prepend,
                self._length_stopper.can_prepend
            ]
        )

        trial_trajectory = partial_trajectory.reversed + shooting_point.trajectory[shooting_point.index + 1:]

        trial_point = paths.ShootingPoint(
            shooting_point.selector,
            trial_trajectory,
            len(partial_trajectory) - 1
        )

        return trial_point

class RandomSubtrajectorySelectGenerator(SampleGenerator):
    '''
    Samples a random subtrajectory satisfying the given subensemble.

    If there are no subtrajectories which satisfy the subensemble, this
    returns the zero-length trajectory.
    '''
    def __init__(self, ensemble, subensemble, n_l):
        super(RandomSubtrajectorySelectGenerator, self).__init__(
            in_ensembles=[ensemble],
            out_ensembles=[subensemble]
        )
        self.n_l = n_l
        self.subensemble = subensemble

    def _choose(self, trajectory_list):
        return random.choice(trajectory_list)

    def __call__(self, trial):
        initial_trajectory = trial.trajectory
        replica = trial.replica
        logger.debug("Working with replica " + str(replica) + " (" + str(initial_trajectory) + ")")

        subtrajs = self.subensemble.split(initial_trajectory)
        logger.debug("Found "+str(len(subtrajs))+" subtrajectories.")

        if (self.n_l is None and len(subtrajs) > 0) or \
            (self.n_l is not None and len(subtrajs) == self.n_l):
            subtraj = self._choose(subtrajs)
        else:
            # return zero-length trajectory otherwise
            subtraj = paths.Trajectory([])

        bias = 1.0

        trial = paths.Sample(
            replica=replica,
            trajectory=subtraj,
            ensemble=self.subensemble,
            parent=trial,
            mover=self,
            bias=bias
        )

        trials = [trial]
        return trials

class ReversalGenerator(SampleGenerator):
    def __call__(self, trial):
        trajectory = trial.trajectory
        ensemble = trial.ensemble
        replica = trial.replica

        reversed_trajectory = trajectory.reversed

        valid = ensemble(reversed_trajectory)
        logger.info("PathReversal move accepted: "+str(valid))

        bias = 1.0

        trial = paths.Sample(
            replica=replica,
            trajectory=reversed_trajectory,
            ensemble=ensemble,
            mover=self,
            parent=trial,
            bias=bias
        )

        return [trial]

class ExtendingGenerator(SampleGenerator):

    def __init__(self, ensemble=None):
        super(ExtendingGenerator, self).__init__()
        self.ensemble = ensemble

    def __call__(self, trial):

        initial_trajectory = trial.trajectory

        if self.ensemble is None:
            dynamics_ensemble = trial.ensemble
        else:
            dynamics_ensemble = self.ensemble

        replica = trial.replica

        trial_trajectory = self._extend(initial_trajectory, dynamics_ensemble)

        trial_details = paths.SampleDetails(
        )

        # the actual bias would be 0.0 since we will never be able to do the
        # reverse move. Since this is the opposite of subtraj we set both
        # proposal bias for these to 100% which means no bias

        trial = paths.Sample(
            replica=replica,
            trajectory=trial_trajectory,
            ensemble=dynamics_ensemble,
            parent=trial,
            details=trial_details,
            mover=self,
            bias=1.0
        )

        trials = [trial]

        return trials

    def _extend(self, initial_trajectory, ensemble):
        return initial_trajectory


class ForwardExtendGenerator(ExtendingGenerator):
    def _shoot(self, initial_trajectory, ensemble):
        shoot_str = "Extending {sh_dir} from frame {fnum} in [0:{maxt}]"
        logger.info(shoot_str.format(fnum=len(initial_trajectory)-1,
                                     maxt=len(initial_trajectory)-1,
                                     sh_dir="forward",
                                    ))

        # Run until one of the stoppers is triggered
        partial_trajectory = SampleGenerator.engine.generate(
            initial_trajectory[-1],
            running = [
                paths.ForwardAppendedTrajectoryEnsemble(
                    ensemble,
                    initial_trajectory[:-1]
                ).can_append,
                self._length_stopper.can_append
            ]
        )

        trial_trajectory = initial_trajectory + partial_trajectory[1:]

        return trial_trajectory


class BackwardExtendGenerator(ExtendingGenerator):
    def _shoot(self, initial_trajectory, ensemble):
        shoot_str = "Extending {sh_dir} from frame {fnum} in [0:{maxt}]"
        logger.info(shoot_str.format(fnum=0,
                                     maxt=len(initial_trajectory)-1,
                                     sh_dir="backward",
                                    ))

        # Run until one of the stoppers is triggered
        partial_trajectory = SampleGenerator.engine.generate(
            initial_trajectory[0].reversed,
            running = [
                paths.BackwardPrependedTrajectoryEnsemble(
                    ensemble,
                    initial_trajectory[1:]
                ).can_prepend,
                self._length_stopper.can_prepend
            ]
        )

        trial_trajectory = partial_trajectory.reversed + initial_trajectory[1:]

        return trial_trajectory
