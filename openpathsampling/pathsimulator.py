import time
import sys
import logging

from openpathsampling.netcdfplus import StorableNamedObject, StorableObject

import openpathsampling as paths
import openpathsampling.tools

from openpathsampling.pathmover import SubPathMover
from ops_logging import initialization_logging
import abc

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


class MCStep(StorableObject):
    """
    A monte-carlo step in the main PathSimulation loop

    It references all objects created and used in a MC step. The used mover,
    and simulator as well as the initial and final sampleset, the step
    number and the generated pathmovechange.

    Attributes
    ----------
    simulation : PathSimulation
        the running pathsimulation responsible for generating the step
    mccycle : int
        the step number counting from the root sampleset
    previous : SampleSet
        the initial (pre) sampleset
    active : SampleSet
        the final (post) sampleset
    change : PathMoveChange
        the pathmovechange describing the transition from pre to post
    """
    def __init__(self,
                 simulation=None,
                 mccycle=-1,
                 previous=None,
                 active=None,
                 change=None
                 ):

        super(MCStep, self).__init__()
        self.simulation = simulation
        self.previous = previous
        self.active = active
        self.change = change
        self.mccycle = mccycle


class PathSimulator(StorableNamedObject):
    __metaclass__ = abc.ABCMeta

    calc_name = "PathSimulator"
    _excluded_attr = ['globalstate', 'step', 'save_frequency']

    def __init__(self, storage):
        super(PathSimulator, self).__init__()
        self.storage = storage
        # self.engine = engine
        self.save_frequency = 1
        self.step = 0
        initialization_logging(
            logger=init_log, obj=self,
            entries=['storage']#, 'engine']
        )

        self.globalstate = None

    # TODO: Remove, is not used
    def set_replicas(self, samples):
        self.globalstate = paths.SampleSet(samples)

    def sync_storage(self):
        """
        Will sync all collective variables and the storage to disk
        """
        if self.storage is not None:
            self.storage.sync_all()

    @abc.abstractmethod
    def run(self, nsteps):
        """
        Run the simulator for a number of steps

        Parameters
        ----------
        nsteps : int
            number of step to be run
        """
        pass

    def save_initial(self):
        """
        Save the initial state as an MCStep to the storage
        """
        mcstep = MCStep(
            simulation=self,
            mccycle=self.step,
            previous=None,
            active=self.globalstate,
            change=paths.EmptyPathMoveChange()
        )

        if self.storage is not None:
            self.storage.steps.save(mcstep)
            self.storage.sync_all()


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
        self._ensemble_dict = {ens : rep for rep, ens in enumerate(ensembles) }
        
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

        self.globalstate = paths.SampleSet([sample])

        if movers is None:
            pass # TODO: implement defaults: one per ensemble, uniform sel
        else:
            self.movers = movers
        initialization_logging(init_log, self,
                               ['movers', 'ensembles'])
        init_log.info("Parameter: %s : %s", 'trajectory', str(trajectory))

        self._bootstrapmove = BootstrapPromotionMove(bias=None,
                                               shooters=self.movers,
                                               ensembles=self.ensembles
                                              )


    def run(self, nsteps):
        bootstrapmove = self._bootstrapmove

        cvs = []
        n_samples = 0

        if self.storage is not None:
            cvs = list(self.storage.cvs)
            n_samples = len(self.storage.snapshots)

        ens_num = len(self.globalstate)-1

        if self.step == 0:
            self.save_initial()

        failsteps = 0
        # if we fail nsteps times in a row, kill the job

        while ens_num < len(self.ensembles) - 1 and failsteps < nsteps:
            self.step += 1
            logger.info("Step: " + str(self.step)
                        + "   Ensemble: " + str(ens_num)
                        + "  failsteps = " + str(failsteps)
                       )
            paths.tools.refresh_output(
                ("Working on Bootstrapping cycle step %d" +
                " in ensemble %d/%d .\n") %
                ( self.step, ens_num + 1, len(self.ensembles) )
            )

            movepath = bootstrapmove.move(self.globalstate)
            samples = movepath.results
            new_sampleset = self.globalstate.apply_samples(samples)

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
                previous=self.globalstate,
                active=new_sampleset,
                change=movepath
            )


#            logger.debug("GLOBALSTATE:")
#            for sample in self.globalstate:
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

            self.globalstate = new_sampleset

            old_ens_num = ens_num
            ens_num = len(self.globalstate)-1
            if ens_num == old_ens_num:
                failsteps += 1

            if self.step % self.save_frequency == 0:
                self.globalstate.sanity_check()
                self.sync_storage()

        self.sync_storage()

        paths.tools.refresh_output(
                ("DONE! Completed Bootstrapping cycle step %d" +
                " in ensemble %d/%d .\n") %
                ( self.step, ens_num + 1, len(self.ensembles) )
            )


class FullBootstrapping(PathSimulator):
    """
    Takes a snapshot as input; gives you back a sample set with trajectories
    for every ensemble in the transition.

    Someday this will be combined with the regular bootstrapping code. 
    """
    calc_name = "FullBootstrapping"

    def __init__(self, transition, snapshot, storage=None, engine=None,
                 extra_interfaces=None, forbidden_states=None, initial_max_length=None):
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
            self.first_traj_ensemble = paths.LengthEnsemble(slice(0, self.initial_max_length)) & self.first_traj_ensemble

        self.extra_ensembles = [paths.TISEnsemble(transition.stateA,
                                                  transition.stateB, iface,
                                                  transition.orderparameter)
                                for iface in extra_interfaces
        ]

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


    def run(self, max_ensemble_rounds=None, n_steps_per_round=20, build_attempts = 20):
        #print first_traj_ensemble #DEBUG
        has_AA_path = False
        while not has_AA_path:
            self.engine.current_snapshot = self.snapshot.copy()
            self.engine.snapshot = self.snapshot.copy()
            print "Building first trajectory"
            sys.stdout.flush()
            first_traj = self.engine.generate(
                self.engine.current_snapshot, 
                [self.first_traj_ensemble.can_append]
            )
            print "Selecting segment"
            sys.stdout.flush()
            subtrajs = self.ensemble0.split(first_traj)
            if len(subtrajs) > 0:
                # if we have a short enough path go ahead
                subtraj = subtrajs[0]
                # check that this is A->A as well
                has_AA_path = self.state(subtraj[-1]) and self.state(subtraj[0])

            build_attempts -= 1
            if build_attempts == 0:
                raise RuntimeError('Too many attempts. Try another initial snapshot instead.')

            
        print "Sampling " + str(self.n_ensembles) + " ensembles."
        bootstrap = paths.Bootstrapping(
            storage=self.storage,
            ensembles=self.all_ensembles,
            movers=self.transition_shooters + self.extra_shooters,
            trajectory=subtraj
        )
        print "Beginning bootstrapping"
        n_rounds = 0
        n_filled = len(bootstrap.globalstate)
        while n_filled < self.n_ensembles:
            bootstrap.run(n_steps_per_round)

            if n_filled == len(bootstrap.globalstate):
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
                else:
                    logger.warning(msg)
                break
            n_filled = len(bootstrap.globalstate)

        return bootstrap.globalstate


class PathSampling(PathSimulator):
    """
    General path sampling code. 
    
    Takes a single move_scheme and generates samples from that, keeping one
    per replica after each move. 
    """

    calc_name = "PathSampling"
    def __init__(
            self,
            storage,
            move_scheme=None,
            globalstate=None
    ):
        """
        Parameters
        ----------
        storage : openpathsampling.storage.Storage
            the storage where all results should be stored in
        engine : openpathsampling.DynamicsEngine
            the engine to be used with shooting moves
        move_scheme : openpathsampling.MoveScheme
            the move scheme used for the pathsampling cycle
        globalstate : openpathsampling.SampleSet
            the initial SampleSet for the Simulator
        """
        super(PathSampling, self).__init__(storage)
        self.move_scheme = move_scheme
        self.root_mover = move_scheme.move_decision_tree()
#        self.move_scheme.name = "PathSamplingRoot"

        samples = []
        if globalstate is not None:
            for sample in globalstate:
                samples.append(sample.copy_reset())

        self.globalstate = paths.SampleSet(samples)
        self.root = self.globalstate

        initialization_logging(init_log, self, 
                               ['move_scheme', 'globalstate'])
        self.live_visualization = None
        self.visualize_frequency = 1
        self._mover = paths.PathSimulatorMover(self.root_mover, self)

    def run_until(self, nsteps):
        if self.storage is not None:
            if len(self.storage.steps) > 0:
                self.step = len(self.storage.steps)
        nsteps_to_run = nsteps - self.step
        self.run(nsteps_to_run)

    def run(self, nsteps):
        mcstep = None

        cvs = list()
        n_samples = 0

        if self.storage is not None:
            n_samples = len(self.storage.snapshots)
            cvs = list(self.storage.cvs)

        if self.step == 0:
            if self.storage is not None:
                self.storage.save(self.move_scheme)
            self.save_initial()

        for nn in range(nsteps):
            self.step += 1
            logger.info("Beginning MC cycle " + str(self.step))
            refresh=True
            if self.step % self.visualize_frequency == 0:
                # do we visualize this step?
                if self.live_visualization is not None and mcstep is not None:
                    # do we visualize at all?
                    self.live_visualization.draw_ipynb(mcstep)
                    refresh=False

                paths.tools.refresh_output(
                    "Working on Monte Carlo cycle number " + str(self.step)
                    + ".\n", 
                    refresh=refresh
                )

            time_start = time.time() 
            movepath = self._mover.move(self.globalstate, step=self.step)
            samples = movepath.results
            new_sampleset = self.globalstate.apply_samples(samples)
            time_elapsed = time.time() - time_start

            # TODO: we can save this with the MC steps for timing? The bit
            # below works, but is only a temporary hack
            setattr(movepath.details, "timing", time_elapsed)

            mcstep = MCStep(
                simulation=self,
                mccycle=self.step,
                previous=self.globalstate,
                active=new_sampleset,
                change=movepath
            )


            if self.storage is not None:
                for cv in cvs:
                    n_len = len(self.storage.snapshots)
                    cv(self.storage.snapshots[n_samples:n_len])
                    n_samples = n_len

                self.storage.steps.save(mcstep)

            if self.step % self.save_frequency == 0:
                self.globalstate.sanity_check()
                self.sync_storage()

            self.globalstate = new_sampleset

        self.sync_storage()

        if self.live_visualization is not None and mcstep is not None:
            self.live_visualization.draw_ipynb(mcstep)
        paths.tools.refresh_output(
            "DONE! Completed " + str(self.step) + " Monte Carlo cycles.\n",
            refresh=False
        )

class CommittorSimulation(PathSimulator):
    def __init__(self, storage, engine=None, states=None, randomizer=None,
                 initial_snapshots=None, direction=None):
        super(CommittorSimulation, self).__init__(storage)
        self.engine = engine
        paths.EngineMover.default_engine = engine
        self.states = states
        self.randomizer = randomizer
        try:
            initial_snapshots = list(initial_snapshots)
        except TypeError:
            initial_snapshots = [initial_snapshots]
        self.initial_snapshots = initial_snapshots
        self.direction = direction

        all_state_volume = paths.join_volumes(states)

        # we should always start from a single frame not in any state
        self.starting_ensemble = (
            paths.AllOutXEnsemble(all_state_volume) &
            paths.LengthEnsemble(1)
        )
        # shoot forward until we hit a state
        self.forward_ensemble = paths.SequentialEnsemble([
            paths.AllOutXEnsemble(all_state_volume),
            paths.AllInXEnsemble(all_state_volume) & paths.LengthEnsemble(1)
        ])
        # or shoot backward until we hit a state
        self.backward_ensemble = paths.SequentialEnsemble([
            paths.AllInXEnsemble(all_state_volume) & paths.LengthEnsemble(1),
            paths.AllOutXEnsemble(all_state_volume)
        ])

        self.forward_mover = paths.ForwardExtendMover(
            ensemble=self.starting_ensemble,
            target_ensemble=self.forward_ensemble
        )
        self.backward_mover = paths.BackwardExtendMover(
            ensemble=self.starting_ensemble,
            target_ensemble=self.backward_ensemble
        )

        if self.direction is None:
            self.mover = paths.RandomChoiceMover([self.forward_mover,
                                                  self.backward_mover])
        elif self.direction > 0:
            self.mover = self.forward_mover
        elif self.direction < 0:
            self.mover = self.backward_mover

    def run(self, n_per_snapshot, as_chain=False):
        self.step = 0
        for snapshot in self.initial_snapshots:
            start_snap = snapshot
            # do what we need to get the snapshot set up
            for step in range(n_per_snapshot):
                if as_chain:
                    start_snap = self.randomizer(start_snap)
                else:
                    start_snap = self.randomizer(snapshot)

                sample_set = paths.SampleSet([
                    paths.Sample(replica=0,
                                 trajectory=paths.Trajectory([start_snap]),
                                 ensemble=self.starting_ensemble)
                ])
                sample_set.sanity_check()
                new_pmc = self.mover.move(sample_set)
                samples = new_pmc.results
                new_sample_set = sample_set.apply_samples(samples)

                mcstep = MCStep(
                    simulation=self,
                    mccycle = self.step,
                    previous=sample_set,
                    active=new_sample_set,
                    change=new_pmc
                )

                if self.storage is not None:
                    self.storage.steps.save(mcstep)
                    if self.step % self.save_frequency == 0:
                        self.sync_storage()

                pass


