from openpathsampling.todict import OPSNamed, OPSObject
import openpathsampling as paths
import openpathsampling.tools
from openpathsampling.pathmover import SubPathMover

import logging
from ops_logging import initialization_logging
logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')

class MCStep(OPSObject):
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

class MCObserver(object):
    def __init__(self):
        self.scope = None

    def start(self):
        pass

    def pre(self):
        pass

    def post(self):
        pass

    def final(self):
        pass

class Reporter(MCObserver):
    def __init__(self, start=None, pre=None, post=None, final=None):
        super(Reporter, self).__init__()
        self.fnc_pre = pre
        self.fnc_final = final
        self.fnc_start = start
        self.fnc_post = post

    def _write(self, s):
        print s

    def start(self):
        self._write(
            self.fnc_start(self.scope)
        )

    def final(self):
        self._write(
            self.fnc_final(self.scope)
        )

    def pre(self):
        self._write(
            self.fnc_pre(self.scope)
        )

    def post(self):
        self._write(
            self.fnc_post(self.scope)
        )


class InfoLoggingReporter(Reporter):
    def _write(self, s):
        logger.info(s)

class IPythonReporter(Reporter):
    def _write(self, s):
        paths.tools.refresh_output(s)


class PathSimulator(OPSNamed):

    _excluded_attr = ['globalstate', 'step', 'save_frequency']

    def __init__(self, storage, engine=None):
        super(PathSimulator, self).__init__()
        self.storage = storage
        self.engine = engine
        self.save_frequency = 1
        self.step = 0
        initialization_logging(
            logger=init_log, obj=self,
            entries=['storage', 'engine']
        )

        self.globalstate = None
        self.last_mcstep = None
        self._mover = None
        self.observer = []

    def __add__(self, other):
        if isinstance(other, MCObserver):
            other.scope = self
            self.observer.append(other)

    def sync_storage(self):
        """
        Will sync all collective variables and the storage to disk
        """
        if self.storage is not None:
            self.storage.cvs.sync()
            self.storage.sync()

    def _is_running(self, n_steps):
        """
        Returns True if the simulation should keep running

        Default is to run until the number of steps is reached. Can be used to
        override the default behaviour and run until a specific event occurs like
        reaching a state, etc...
        """
        return self.current_step < n_steps

    def run(self, n_steps):
        """
        Run the simulator for a number of steps

        Parameters
        ----------
        nsteps : int
            number of step to be run
        """
        mcstep = None

        cvs = list()
        n_samples = 0

        if self.storage is not None:
            n_samples = len(self.storage.snapshots)
            cvs = list(self.storage.cvs)

        if self.step == 0:
            self.save_initial()

        self.current_step = 0

        for obs in self.observer.values():
            obs.start()

        while self._is_running(n_steps):
            self.step += 1
            self.current_step += 1

            for obs in self.observer.values():
                obs.pre()

            movepath = self._mover.move(self.globalstate, step=self.step)
            samples = movepath.results
            new_sampleset = self.globalstate.apply_samples(samples)

            self.last_mcstep = MCStep(
                simulation=self,
                mccycle=self.step,
                previous=self.globalstate,
                active=new_sampleset,
                change=movepath
            )

            for obs in self.observer.values():
                obs.post()

            if self.storage is not None:
                for cv in cvs:
                    n_len = len(self.storage.snapshots)
                    cv(self.storage.snapshots[n_samples:n_len])
                    n_samples = n_len

                self.storage.steps.save(self.last_mcstep)

            if self.step % self.save_frequency == 0:
                self.globalstate.sanity_check()
                self.sync_storage()

            self.globalstate = new_sampleset

        self.sync_storage()

        for obs in self.observer.values():
            obs.final()

    def save_initial(self):
        """
        Save the initial state as an MCStep to the storage
        """
        self.last_mcstep = MCStep(
            simulation=self,
            mccycle=self.step,
            previous=None,
            active=self.globalstate,
            change=paths.EmptyPathMoveChange()
        )

        if self.storage is not None:
            self.storage.steps.save(self.last_mcstep)
            self.storage.sync()


class BootstrapPromotionMove(SubPathMover):
    """
    Bootstrap promotion is the combination of an EnsembleHop (to the next
    ensemble up) with incrementing the replica ID.
    """
    def __init__(self, bias=None, shooters=None,
                 ensembles=None):
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
        movers : list of openpathsampling.PathMover or None
            list of shooters to be used in the BootstrapPromotionMove. If None then
            Bootstrapping creates uniform selecting OneWayShooters for each ensemble.
        trajectory : openpathsampling.Trajectory
            an initial trajectory to be started from
        ensembles : nested list of openpathsampling.Ensemble
            the ensembles this move should act on
        """
        # TODO: Change input from trajectory to sample
        super(Bootstrapping, self).__init__(storage, engine)
        self.ensembles = ensembles

        sample = paths.Sample(
            replica=0,
            trajectory=trajectory,
            ensemble=self.ensembles[0]
        )

        self.globalstate = paths.SampleSet([sample])
        if self.storage is not None:
            self.storage.samplesets.save(self.globalstate)
        if movers is None:
            movers = [
                paths.OneWayShootingMover(
                    ensemble=ens,
                    selector=paths.UniformSelector()
                )
                for ens in self.ensembles
            ]
        else:
            self.movers = movers

        initialization_logging(init_log, self,
                               ['movers', 'ensembles'])
        init_log.info("Parameter: %s : %s", 'trajectory', str(trajectory))

        self._mover = BootstrapPromotionMove(bias=None,
                                               shooters=self.movers,
                                               ensembles=self.ensembles
                                              )

        self.failsteps = 0
        self.old_ens_num = 0

        self + IPythonReporter(
            pre = lambda this :
                ("Working on Bootstrapping cycle step %d" +
                    " in ensemble %d/%d .\n") %
                    ( this.step, len(this.globalstate), len(this.ensembles) ),
            final = lambda this :
                ("DONE! Completed Bootstrapping cycle step %d" +
                    " in ensemble %d/%d .\n") %
                    ( this.step, len(this.globalstate), len(this.ensembles) )
        )

        self + InfoLoggingReporter(
            pre = lambda this : "Beginning Bootstrapping cycle " + str(this.step)
        )

    def _is_running(self, n_steps):

        ens_num = len(self.globalstate) - 1
        if ens_num == self.old_ens_num:
            self.failsteps += 1
        else:
            self.old_ens_num = ens_num

        return len(self.globalstate) < len(self.ensembles) and self.failsteps < n_steps


class PathSampling(PathSimulator):
    """
    General path sampling code. 
    
    Takes a single move_scheme and generates samples from that, keeping one
    per replica after each move. 
    """

    def __init__(
            self,
            storage,
            engine=None,
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
        move_scheme : openpathsampling.PathMover
            the mover used for the pathsampling cycle
        globalstate : openpathsampling.SampleSet
            the initial SampleSet for the Simulator
        """
        super(PathSampling, self).__init__(storage, engine)
        self.move_scheme = move_scheme
#        self.move_scheme.name = "PathSamplingRoot"

        samples = []
        if globalstate is not None:
            for sample in globalstate:
                samples.append(sample.copy_reset())

        self.globalstate = paths.SampleSet(samples)

        initialization_logging(init_log, self, 
                               ['move_scheme', 'globalstate'])

        self._mover = paths.PathSimulatorMover(self.move_scheme, self)

        self + IPythonReporter(
            pre = lambda this :  "Working on Monte Carlo cycle step " + str(this.step) + ".\n",
            final = lambda this : "DONE! Completed " + str(this.step) + " Monte Carlo cycles.\n"
        )

        self + InfoLoggingReporter(
            pre = lambda this : "Beginning MC cycle " + str(this.step)
        )