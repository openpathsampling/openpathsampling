from pathmover import (PathMover, MoveDetails, ReplicaExchangeMover,
                       EnsembleHopMover)
from sample import SampleSet, Sample

from openpathsampling.todict import restores_as_stub_object
from openpathsampling.pathmover import PathMover

from openpathsampling.movepath import SampleMovePath

import openpathsampling as paths

import logging
from ops_logging import initialization_logging
logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')

class Calculation(object):

    calc_name = "Calculation"

    def __init__(self, storage, engine=None):
        self.storage = storage
        self.engine = engine
        initialization_logging(
            logger=init_log, obj=self,
            entries=['storage', 'engine']
        )

    def set_replicas(self, samples):
        self.globalstate = SampleSet(samples)

    def run(self, nsteps):
        logger.warning("Running an empty calculation? Try a subclass, maybe!")


@restores_as_stub_object
class BootstrapPromotionMove(PathMover):
    '''
    Bootstrap promotion is the combination of an EnsembleHop (to the next
    ensemble up) with incrementing the replica ID.
    '''
    def __init__(self, bias=None, shooters=None,
                 ensembles=None, replicas='all'):
        super(BootstrapPromotionMove, self).__init__(ensembles=ensembles,
                                                     replicas=replicas)
        self.shooters = shooters
        self.bias = bias
        initialization_logging(logger=init_log, obj=self,
                               entries=['bias', 'shooters'])


        # Create all possible hoppers so we do not have to recreate these
        # every time which will result in more efficient storage

        ens_pairs = []

        for rep in range(0,len(self.ensembles) - 1):
            ensemble_from = self.ensembles[rep]
            ensemble_to = self.ensembles[rep+1]

            ens_pairs.append([ensemble_from, ensemble_to])

        self._hopper = EnsembleHopMover(
            bias=self.bias,
            ensembles=ens_pairs,
            replicas='all'
        )


    def move(self, globalstate):
        # the tricky part here is that, if the hop is allowed, we only want
        # to report the sample in the new ensemble. The way we do this is by
        # treating each bootstrap move as a combination of 3 moves: it
        # always starts with a shooting move and a replica hop, and then, if
        # the hop was successful, a replica ID change move

        top_rep = max(globalstate.replica_list())

        shooter = self.shooters[top_rep]

        mover = paths.PartialAcceptanceSequentialMover([
            shooter,
            self._hopper
        ])

        return mover.move(globalstate)

class InitializeSingleTrajectoryMover(PathMover):
    def __init__(self, bias=None, shooters=None,
                 ensembles=None, replicas='all'):
        super(InitializeSingleTrajectoryMover, self).__init__(ensembles=ensembles,
                                                     replicas=replicas)
        self.shooters = shooters
        self.bias = bias
        initialization_logging(logger=init_log, obj=self,
                               entries=['bias', 'shooters'])

    def move(self, globalstate=None):
        init_details = MoveDetails()
        init_details.accepted = True
        init_details.acceptance_probability = 1.0
        init_details.mover = self
        init_details.inputs = []
        init_details.trial = trajectory
        init_details.ensemble = ensemble
        sample = Sample(replica=0, trajectory=trajectory,
                        ensemble=self.ensembles[0], details=init_details)


class Bootstrapping(Calculation):
    """The ensembles for the Bootstrapping calculation must be one ensemble
    set, in increasing order."""

    calc_name = "Bootstrapping"

    def __init__(self, storage, engine=None, movers=None, trajectory=None,
                 ensembles=None):
        super(Bootstrapping, self).__init__(storage, engine)
        self.ensembles = ensembles

        # this is stupid; must be a better way
        init_details = MoveDetails()
        init_details.accepted = True
        init_details.acceptance_probability = 1.0
        init_details.mover = PathMover()
        init_details.mover.name = "Initialization (trajectory)"
        init_details.inputs = []
        init_details.trial = trajectory
        init_details.ensemble = self.ensembles[0]
        sample = Sample(replica=0, trajectory=trajectory, 
                        ensemble=self.ensembles[0], details=init_details)

        self.globalstate = SampleSet([sample])
        if self.storage is not None:
            self.globalstate.save_samples(self.storage)
        if movers is None:
            pass # TODO: implement defaults: one per ensemble, uniform sel
        else:
            self.movers = movers
        initialization_logging(init_log, self,
                               ['movers', 'ensembles'])
        init_log.info("Parameter: %s : %s", 'trajectory', str(trajectory))

    def run(self, nsteps):
        # TODO: turn off init_log during run loop
        bootstrapmove = BootstrapPromotionMove(bias=None,
                                               shooters=self.movers,
                                               ensembles=self.ensembles,
                                               replicas='all'
                                              )

        ens_num = 0
        failsteps = 0
        step_num = 0
        # if we fail nsteps times in a row, kill the job

        old_ens = self.globalstate[0].ensemble

        while ens_num < len(self.ensembles) - 1 and failsteps < nsteps:
            logger.info("Step: " + str(step_num) 
                        + "   Ensemble: " + str(ens_num)
                        + "  failsteps = " + str(failsteps)
                       )

            movepath = bootstrapmove.move(self.globalstate)
            samples = movepath.samples
            self.globalstate = self.globalstate.apply_samples(samples, step=step_num)

            if movepath.movepaths[0].accepted is True:
                # shooter has been accepted
                if movepath.movepaths[1].accepted is True:
                    # hop has been accepted!
                    failsteps = 0
                else:
                    failsteps += 1

            if self.storage is not None:
                self.globalstate.save_samples(self.storage)
                self.globalstate.save(self.storage)
                # TODO NEXT: store self.globalstate itself
            step_num += 1

        for sample in self.globalstate:
            assert sample.ensemble(sample.trajectory) == True, "WTF?"

class PathSampling(Calculation):
    """
    General path sampling code. Takes a single root_mover and generates
    samples from that, keeping one per replica after each move. 

    TODO
    ----
        Add a nice syntax for the root_mover, at least allowing us to
        implicitly combine simultaneous multimovers with mixed movers.
    """

    calc_name = "PathSampling"
    def __init__(self, storage, engine=None, root_mover=None,
                 globalstate=None):
        super(PathSampling, self).__init__(storage, engine)
        self.root_mover = root_mover
        self.root_mover.name = "PathSamplingRoot"
        samples = []
        for sample in globalstate:
            samples.append(sample.copy_reset())
        self.globalstate = SampleSet(samples)

        initialization_logging(init_log, self, 
                               ['root_mover', 'globalstate'])


    def run(self, nsteps):
        # TODO: change so we can start from some arbitrary step number
    
        self.globalstate = self.globalstate.apply_samples(self.globalstate,
                                                          step=-1)

        if self.storage is not None:
            self.globalstate.save_samples(self.storage)
            self.globalstate.save(self.storage)
            self.storage.sync()

        for step in range(nsteps):
            samples = self.root_mover.move(self.globalstate)
            self.globalstate = self.globalstate.apply(samples,
                                                              step=step)
            if self.storage is not None:
                self.globalstate.save_samples(self.storage)
                self.globalstate.save(self.storage)
                self.storage.sync()
