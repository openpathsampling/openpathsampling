import sys
import collections
import logging

import openpathsampling as paths

from openpathsampling.pathmover import SubPathMover
from .path_simulator import PathSimulator, MCStep
from ..ops_logging import initialization_logging

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')

class BootstrapPromotionMove(SubPathMover):
    """
    Bootstrap promotion is the combination of an EnsembleHop (to the next
    ensemble up) with incrementing the replica ID.
    """
    def __init__(self, bias=None, shooters=None, ensembles=None):
        """
        Parameters
        ----------
        bias : None
            not used yet, only for API consistency and later implementation
        shooters : list of ShootingMovers
            list of ShootingMovers for each ensemble
        ensembles : list of Ensembles
            list of ensembles the move should act on

        Notes
        -----
        The bootstrapping will use the ensembles sequentially so it requires
        that all ensembles have a reasonable overlab using shooting moves.

        """
        self.shooters = shooters
        self.bias = bias
        self.ensembles = ensembles
        initialization_logging(logger=init_log, obj=self,
                               entries=['bias', 'shooters', 'ensembles'])

        ens_pairs = [[self.ensembles[i], self.ensembles[i+1]]
                     for i in range(len(self.ensembles)-1)]

        # Bootstrapping sets numeric replica IDs. If the user wants it done
        # differently, the user can change it.
        self._ensemble_dict = {ens: rep for rep, ens in enumerate(ensembles)}

        # Create all possible hoppers so we do not have to recreate these
        # every time which will result in more efficient storage
        mover = paths.LastAllowedMover([
            # writing an algorithm this convoluted can get you shot in Texas
            paths.PartialAcceptanceSequentialMover(
                movers=[
                    shoot,
                    paths.EnsembleHopMover(
                        ensemble=enss[0],
                        target_ensemble=enss[1],
                        change_replica=self._ensemble_dict[enss[1]]
                    )
                ]
            ) for (enss, shoot) in zip(ens_pairs, shooters)
        ])

        super(BootstrapPromotionMove, self).__init__(mover)


class Bootstrapping(PathSimulator):
    """Creates a SampleSet with one sample per ensemble.

    The ensembles for the Bootstrapping pathsimulator must be one ensemble
    set, in increasing order. Replicas are named numerically.
    """

    calc_name = "Bootstrapping"

    def __init__(
            self,
            storage,
            engine=None,
            movers=None,
            trajectory=None,
            ensembles=None
    ):
        """
        Parameters
        ----------
        storage : openpathsampling.storage.Storage
            the storage all results should be stored in
        engine : openpathsampling.DynamicsEngine
            the dynamics engine to be used
        movers : list of openpathsampling.PathMover
            list of shooters to be used in the BootstrapPromotionMove
        trajectory : openpathsampling.Trajectory
            an initial trajectory to be started from
        ensembles : nested list of openpathsampling.Ensemble
            the ensembles this move should act on
        """
        # TODO: Change input from trajectory to sample
        super(Bootstrapping, self).__init__(storage)
        self.engine = engine
        paths.EngineMover.default_engine = engine  # set the default
        self.ensembles = ensembles
        self.trajectory = trajectory

        sample = paths.Sample(
            replica=0,
            trajectory=trajectory,
            ensemble=self.ensembles[0]
        )

        self.sample_set = paths.SampleSet([sample])

        if movers is None:
            pass # TODO: implement defaults: one per ensemble, uniform sel
        else:
            self.movers = movers
        initialization_logging(init_log, self,
                               ['movers', 'ensembles'])
        init_log.info("Parameter: %s : %s", 'trajectory', str(trajectory))

        self._bootstrapmove = BootstrapPromotionMove(
            bias=None,
            shooters=self.movers,
            ensembles=self.ensembles
        )

    def run(self, n_steps):
        bootstrapmove = self._bootstrapmove

        cvs = []
        n_samples = 0

        if self.storage is not None:
            cvs = list(self.storage.cvs)
            n_samples = len(self.storage.snapshots)

        ens_num = len(self.sample_set)-1

        if self.step == 0:
            self.save_initial_step()

        failsteps = 0
        # if we fail n_steps times in a row, kill the job

        while ens_num < len(self.ensembles) - 1 and failsteps < n_steps:
            self.step += 1
            logger.info("Step: " + str(self.step)
                        + "   Ensemble: " + str(ens_num)
                        + "  failsteps = " + str(failsteps)
                       )
            paths.tools.refresh_output(
                ("Working on Bootstrapping cycle step %d" +
                 " in ensemble %d/%d .\n") %
                (self.step, ens_num + 1, len(self.ensembles)),
                output_stream=self.output_stream,
                refresh=self.allow_refresh
            )

            movepath = bootstrapmove.move(self.sample_set)
            samples = movepath.results
            new_sampleset = self.sample_set.apply_samples(samples)

#            samples = movepath.results
#            logger.debug("SAMPLES:")
#            for sample in samples:
#                logger.debug("(" + str(sample.replica)
#                             + "," + str(sample.trajectory)
#                             + "," + repr(sample.ensemble)
#                            )


            mcstep = MCStep(
                simulation=self,
                mccycle=self.step,
                previous=self.sample_set,
                active=new_sampleset,
                change=movepath
            )


#            logger.debug("GLOBALSTATE:")
#            for sample in self.sample_set:
#                logger.debug("(" + str(sample.replica)
#                             + "," + str(sample.trajectory)
#                             + "," + repr(sample.ensemble)
#                            )



            if self.storage is not None:
                # compute all cvs now
                for cv in cvs:
                    n_len = len(self.storage.snapshots)
                    cv(self.storage.snapshots[n_samples:n_len])
                    n_samples = n_len

                self.storage.steps.save(mcstep)

            self.sample_set = new_sampleset

            old_ens_num = ens_num
            ens_num = len(self.sample_set)-1
            if ens_num == old_ens_num:
                failsteps += 1

            if self.step % self.save_frequency == 0:
                self.sample_set.sanity_check()
                self.sync_storage()

        self.sync_storage()

        paths.tools.refresh_output(
            ("DONE! Completed Bootstrapping cycle step %d"
             + " in ensemble %d/%d.\n") %
            (self.step, ens_num + 1, len(self.ensembles)),
            output_stream=self.output_stream,
            refresh=self.allow_refresh
        )


class FullBootstrapping(PathSimulator):
    """
    Takes a snapshot as input; gives you back a sampleset with trajectories
    for every ensemble in the transition.

    This includes

    Parameters
    ----------
    transition : :class:`.TISTransition`
        the TIS transition to fill by bootstrapping
    snapshot : :class:`.Snapshot`
        the initial snapshot
    storage : :class:`.Storage`
        storage file to record the steps (optional)
    engine : :class:`.DynamicsEngine`
        MD engine to use for dynamics
    extra_interfaces : list of :class:`.Volume`
        additional interfaces to make into TIS ensembles (beyond those in
        the transition)
    extra_ensembles : list of :class:`.Ensemble`
        additional ensembles to sample after the TIS ensembles
    forbidden_states : list of :class:`.Volume`
        regions that are disallowed during the initial trajectory. Note that
        these region *are* allowed during the interface sampling
    initial_max_length : int
        maximum length of the initial A->A trajectory
    """
    calc_name = "FullBootstrapping"

    def __init__(self, transition, snapshot, storage=None, engine=None,
                 extra_interfaces=None, extra_ensembles=None,
                 forbidden_states=None, initial_max_length=None):
        super(FullBootstrapping, self).__init__(storage)
        self.engine = engine
        paths.EngineMover.default_engine = engine  # set the default
        if extra_interfaces is None:
            extra_interfaces = list()

        if forbidden_states is None:
            forbidden_states = list()
        interface0 = transition.interfaces[0]
        ensemble0 = transition.ensembles[0]
        state = transition.stateA
        self.state = state
        self.first_traj_ensemble = paths.SequentialEnsemble([
            paths.OptionalEnsemble(paths.AllOutXEnsemble(state)),
            paths.AllInXEnsemble(state),
            paths.OptionalEnsemble(
                paths.AllOutXEnsemble(state) & paths.AllInXEnsemble(interface0)
            ),
            paths.OptionalEnsemble(paths.AllInXEnsemble(interface0)),
            paths.AllOutXEnsemble(interface0),
            paths.OptionalEnsemble(paths.AllOutXEnsemble(state)),
            paths.SingleFrameEnsemble(paths.AllInXEnsemble(state))
        ]) & paths.AllOutXEnsemble(paths.join_volumes(forbidden_states))

        self.initial_max_length = initial_max_length

        if self.initial_max_length is not None:
            self.first_traj_ensemble = (
                paths.LengthEnsemble(slice(0, self.initial_max_length))
                & self.first_traj_ensemble
            )

        if extra_ensembles is None:
            extra_ensembles = []
        self.extra_ensembles = [
            paths.TISEnsemble(transition.stateA, transition.stateB, iface,
                              transition.orderparameter)
            for iface in extra_interfaces
        ] + extra_ensembles

        self.transition_shooters = [
            paths.OneWayShootingMover(selector=paths.UniformSelector(),
                                      ensemble=ens,
                                      engine=self.engine)
            for ens in transition.ensembles
        ]

        self.extra_shooters = [
            paths.OneWayShootingMover(selector=paths.UniformSelector(),
                                      ensemble=ens,
                                      engine=self.engine)
            for ens in self.extra_ensembles
        ]
        self.snapshot = snapshot.copy()
        self.ensemble0 = ensemble0
        self.all_ensembles = transition.ensembles + self.extra_ensembles
        self.n_ensembles = len(self.all_ensembles)
        self.error_max_rounds = True


    def run(self, max_ensemble_rounds=None, n_steps_per_round=20,
            build_attempts=20):
        #print first_traj_ensemble #DEBUG
        has_AA_path = False
        subtraj = None
        while not has_AA_path:
            self.engine.current_snapshot = self.snapshot.copy()
            self.engine.snapshot = self.snapshot.copy()
            self.output_stream.write("Building first trajectory\n")
            sys.stdout.flush()
            first_traj = self.engine.generate(
                self.engine.current_snapshot,
                [self.first_traj_ensemble.can_append]
            )
            self.output_stream.write("Selecting segment\n")
            sys.stdout.flush()
            subtrajs = self.ensemble0.split(first_traj)
            if len(subtrajs) > 0:
                # if we have a short enough path go ahead
                subtraj = subtrajs[0]
                # check that this is A->A as well
                has_AA_path = self.state(subtraj[-1]) \
                    and self.state(subtraj[0])

            build_attempts -= 1
            if build_attempts == 0:
                raise RuntimeError(
                    'Too many attempts. Try another initial snapshot instead.')

        self.output_stream.write("Sampling " + str(self.n_ensembles) +
                                 " ensembles.\n")
        bootstrap = paths.Bootstrapping(
            storage=self.storage,
            ensembles=self.all_ensembles,
            movers=self.transition_shooters + self.extra_shooters,
            trajectory=subtraj
        )
        bootstrap.output_stream = self.output_stream
        self.output_stream.write("Beginning bootstrapping\n")
        n_rounds = 0
        n_filled = len(bootstrap.sample_set)
        while n_filled < self.n_ensembles:
            bootstrap.run(n_steps_per_round)

            if n_filled == len(bootstrap.sample_set):
                n_rounds += 1
            else:
                n_rounds = 0
            if n_rounds == max_ensemble_rounds:
                # hard equality instead of inequality so that None gives us
                # effectively infinite (rounds add one at a time
                msg = ("Too many rounds of bootstrapping: " + str(n_rounds)
                       + " round of " + str(n_steps_per_round) + " steps.")
                if self.error_max_rounds:
                    raise RuntimeError(msg)
                else:  # pragma: no cover
                    logger.warning(msg)
                    break
            n_filled = len(bootstrap.sample_set)

        return bootstrap.sample_set
