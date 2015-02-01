from pathmover import (PathMover, MoveDetails, ReplicaExchangeMover,
                       EnsembleHopMover)
from sample import SampleSet, Sample

from opentis.todict import restores_as_stub_object
from opentis.pathmover import PathMover

import logging
from ops_logging import initialization_logging
logger = logging.getLogger(__name__)
init_log = logging.getLogger('opentis.initialization')

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

    def move(self, globalstate):
        # the tricky part here is that, if the hop is allowed, we only want
        # to report the sample in the new ensemble. The way we do this is by
        # treating each bootstrap move as a combination of 3 moves: it
        # always starts with a shooting move and a replica hop, and then, if
        # the hop was successful, a replica ID change move

        #print "Starting BootstrapPromotionMove"
        
        # We make extra variables here so that we can easily refactor. The
        # speed cost the negligible, and it makes it easy to change things.
        top_rep = max(globalstate.replica_list())
        ensemble_from = self.ensembles[top_rep]
        ensemble_to = self.ensembles[top_rep+1]
        old_sample = globalstate[top_rep]

        details = MoveDetails()
        details_inputs = [old_sample.trajectory]
        details_mover = self
        details_replica = top_rep

        init_sample_set = SampleSet([old_sample])

        shooter = self.shooters[top_rep]
        # TODO: should move this to normal initialization so we don't init
        # it every time
        hopper = EnsembleHopMover(bias=self.bias,
                                  ensembles=[ensemble_from, ensemble_to],
                                  replicas=top_rep)

        shoot_samp = shooter.move(init_sample_set)[0]
        init_sample_set = init_sample_set.apply_samples(shoot_samp)
        hop_samp = hopper.move(init_sample_set)[0]
        init_sample_set = init_sample_set.apply_samples(hop_samp)

        # bring all the metadata from the submoves into our details
        details.__dict__.update(shoot_samp.details.__dict__)
        details.__dict__.update(hop_samp.details.__dict__)

        # set the rest of the details to their correct values
        details.inputs = details_inputs
        details.mover = details_mover
        details.replica = details_replica
        details.trial = shoot_samp.details.trial # TODO: is it, though?
        # what about hop trial when that happens? may be cleanest if we make
        # this truly sequential

        # the move will be accepted if the shooting move is accepted, no
        # matter what
        details.accepted = (shoot_samp.details.accepted or 
                            hop_samp.details.accepted)

        # result trajectory is whatever came out of hop_samp
        details.result = hop_samp.details.result
        #details.result_ensemble = hop_samp.details.result_ensemble

        setattr(details, 'start_replica', details.replica)
        if hop_samp.details.accepted == True:
            setattr(details, 'result_replica', details.replica+1)
        else:
            setattr(details, 'result_replica', details.replica)
        sample = Sample(replica=details.result_replica,
                        ensemble=details.result_ensemble,
                        trajectory=details.result,
                        details=details)

        logger.debug("BootstrapMover: accepted = " + str(details.accepted))
        logger.debug(" Shooting part: accepted = " + str(shoot_samp.details.accepted))
        logger.debug("  Hopping part: accepted = " + str(hop_samp.details.accepted))

        return [sample]


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
        init_details.inputs = [trajectory]
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

        while ens_num < len(self.ensembles) - 1 and failsteps < nsteps:
            logger.info("Step: " + str(step_num) 
                        + "   Ensemble: " + str(ens_num)
                        + "  failsteps = " + str(failsteps)
                       )
            old_rep = max(self.globalstate.replica_list())
            samples = bootstrapmove.move(self.globalstate)
            self.globalstate = self.globalstate.apply_samples(samples, step=step_num)
            #print self.globalstate.samples[0]

            if samples[0].replica == old_rep:
                failsteps += 1
            else:
                failsteps = 0
                ens_num += 1

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
        self.globalstate = globalstate
        initialization_logging(init_log, self, 
                               ['root_mover', 'globalstate'])


    def run(self, nsteps):
        for step in range(nsteps):
            samples = self.root_mover.move(self.globalstate)
            self.globalstate = self.globalstate.apply_samples(samples,
                                                              step=step)
            if self.storage is not None:
                self.globalstate.save_samples(self.storage)
                self.globalstate.save(self.storage)

