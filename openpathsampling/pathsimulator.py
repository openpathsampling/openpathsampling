from openpathsampling.todict import OPSNamed, OPSObject
import openpathsampling as paths
import openpathsampling.tools

from openpathsampling.pathmover import PathMover

import logging
from ops_logging import initialization_logging
logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')

import sys

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

class PathSimulator(OPSNamed):

    calc_name = "PathSimulator"
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

    # TODO: Remove, is not used
    def set_replicas(self, samples):
        self.globalstate = paths.SampleSet(samples)

    def sync_storage(self):
        """
        Will sync all collective variables and the storage to disk
        """
        if self.storage is not None:
            self.storage.cvs.sync()
            self.storage.sync()

    def run(self, nsteps):
        """
        Run the simulator for a number of steps

        Parameters
        ----------
        nsteps : int
            number of step to be run
        """
        logger.warning("Running an empty pathsimulator? Try a subclass, maybe!")

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
            self.storage.sync()

class BootstrapPromotionMove(PathMover):
    '''
    Bootstrap promotion is the combination of an EnsembleHop (to the next
    ensemble up) with incrementing the replica ID.
    '''
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
        super(BootstrapPromotionMove, self).__init__(ensembles=ensembles)
        self.shooters = shooters
        self.bias = bias
        initialization_logging(logger=init_log, obj=self,
                               entries=['bias', 'shooters'])

        ens_pairs = [[self.ensembles[i], self.ensembles[i+1]]
                     for i in range(len(self.ensembles)-1)]

        # Bootstrapping sets numeric replica IDs. If the user wants it done
        # differently, the user can change it.
        self._ensemble_dict = {ens : rep for rep, ens in enumerate(ensembles) }
        
        # Create all possible hoppers so we do not have to recreate these
        # every time which will result in more efficient storage
        self._hopper = {}
        for (enss, shoot) in zip(ens_pairs, shooters):
            rep_from = self._ensemble_dict[enss[0]]
            rep_to = self._ensemble_dict[enss[1]]
            # writing an algorithm this convoluted can get you shot in Texas
            self._hopper[rep_from] = paths.RestrictToLastSampleMover(
                paths.PartialAcceptanceSequentialMover(
                    movers=[
                        shoot,
                        paths.EnsembleHopMover(ensembles=enss),
                        paths.ReplicaIDChangeMover(
                            replica_pairs=[rep_from, rep_to]
                        )
                    ]
                )
            )


    def move(self, globalstate):
        # find latest ensemble in the list
        top_ens_idx = len(globalstate)-1
        mover = self._hopper[top_ens_idx]
        return mover.move(globalstate)



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

        self.root = self.globalstate

    def run(self, nsteps):
        bootstrapmove = self._bootstrapmove

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
        self.root = self.globalstate

        initialization_logging(init_log, self, 
                               ['move_scheme', 'globalstate'])
        self.live_visualization = None
        self._mover = paths.PathSimulatorMover(self.move_scheme, self)

    def run(self, nsteps):
        mcstep = None

        if self.storage is not None:
            n_samples = len(self.storage.snapshots)
            cvs = list(self.storage.cvs)

        if self.step == 0:
            self.save_initial()

        for nn in range(nsteps):
            self.step += 1
            logger.info("Beginning MC cycle " + str(self.step))
            refresh=True
            if self.live_visualization is not None and mcstep is not None:
                self.live_visualization.draw(mcstep)
                refresh=False

            paths.tools.refresh_output(
                "Working on Monte Carlo cycle step " + str(self.step) + ".\n",
                refresh=refresh
            )

            movepath = self._mover.move(self.globalstate, step=self.step)
            samples = movepath.results
            new_sampleset = self.globalstate.apply_samples(samples)

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
        paths.tools.refresh_output(
            "DONE! Completed " + str(self.step) + " Monte Carlo cycles.\n"
        )
