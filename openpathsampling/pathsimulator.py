from openpathsampling.todict import OPSNamed, OPSObject
import openpathsampling as paths
import openpathsampling.tools

from openpathsampling.pathmover import PathMover

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
    step : int
        the step number counting from the root sampleset
    pre : SampleSet
        the initial (pre) sampleset
    post : SampleSet
        the final (post) sampleset
    change : PathMoveChange
        the pathmovechange describing the transition from pre to post
    root : SampleSet
        the first sampleset at step 0.
    """
    def __init__(self,
                 simulation=None,
                 step=-1,
                 previous=None,
                 active=None,
                 change=None
                 ):

        super(MCStep, self).__init__()
        self.simulation = simulation
        self.previous = previous
        self.active = active
        self.change = change
        self.step = step

class PathSimulator(OPSNamed):

    calc_name = "PathSimulator"
    _excluded_attr = ['globalstate']

    def __init__(self, storage, engine=None):
        super(PathSimulator, self).__init__()
        self.storage = storage
        self.engine = engine
        self.save_frequency = 1
        initialization_logging(
            logger=init_log, obj=self,
            entries=['storage', 'engine']
        )

    def set_replicas(self, samples):
        self.globalstate = paths.SampleSet(samples)

    def sync_storage(self):
        if self.storage is not None:
            self.storage.cvs.sync()
            self.storage.sync()

    def run(self, nsteps):
        logger.warning("Running an empty pathsimulator? Try a subclass, maybe!")


class BootstrapPromotionMove(PathMover):
    '''
    Bootstrap promotion is the combination of an EnsembleHop (to the next
    ensemble up) with incrementing the replica ID.
    '''
    def __init__(self, bias=None, shooters=None,
                 ensembles=None):
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

        self.step = 0

        mcstep = MCStep(
            simulation=self,
            step=self.step,
            previous=None,
            active=self.globalstate,
            change=paths.EmptyPathMoveChange()
        )

        if self.storage is not None:
            self.storage.steps.save(mcstep)
            self.storage.sync()


    def run(self, nsteps):
        bootstrapmove = self._bootstrapmove

        ens_num = len(self.globalstate)-1
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
                step=self.step,
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
                self.storage.steps.save(mcstep)

            self.globalstate = new_sampleset

            old_ens_num = ens_num
            ens_num = len(self.globalstate)-1
            if ens_num == old_ens_num:
                failsteps += 1

            if self.step % self.save_frequency == 0:
                self.globalstate.sanity_check()
                self.sync_storage()


class PathSampling(PathSimulator):
    """
    General path sampling code. 
    
    Takes a single root_mover and generates samples from that, keeping one
    per replica after each move. 
    """

    calc_name = "PathSampling"
    def __init__(
            self,
            storage,
            engine=None,
            root_mover=None,
            globalstate=None
    ):
        super(PathSampling, self).__init__(storage, engine)
        self.root_mover = root_mover
        self.root_mover.name = "PathSamplingRoot"

        samples = []
        if globalstate is not None:
            for sample in globalstate:
                samples.append(sample.copy_reset())

        self.globalstate = paths.SampleSet(samples)
        self.root = self.globalstate

        initialization_logging(init_log, self, 
                               ['root_mover', 'globalstate'])

        self._mover = paths.PathSimulatorMover(self.root_mover, self)

        self.step = 0

        mcstep = MCStep(
            simulation=self,
            step=self.step,
            previous=None,
            active=self.globalstate,
            change=paths.EmptyPathMoveChange()
        )

        if self.storage is not None:
            self.storage.steps.save(mcstep)
            self.storage.sync()

    def run(self, nsteps):
        # TODO: change so we can start from some arbitrary step number

        mcstep = None

        for nn in range(nsteps):
            self.step += 1
            logger.info("Beginning MC cycle " + str(self.step))
            paths.tools.refresh_output(
                "Working on Monte Carlo cycle step " + str(self.step) + ".\n"
            )
            movepath = self._mover.move(self.globalstate, step=self.step)
            samples = movepath.results
            new_sampleset = self.globalstate.apply_samples(samples)

            mcstep = MCStep(
                simulation=self,
                step=self.step,
                previous=self.globalstate,
                active=new_sampleset,
                change=movepath
            )

            if self.storage is not None:
                self.storage.steps.save(mcstep)

            if self.step % self.save_frequency == 0:
                self.globalstate.sanity_check()
                self.sync_storage()

            self.globalstate = new_sampleset


        self.sync_storage()
        paths.tools.refresh_output(
            "DONE! Completed " + str(nsteps+1) + " Monte Carlo cycles.\n"
        )

        return mcstep