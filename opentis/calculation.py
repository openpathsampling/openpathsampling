from globalstate import GlobalState
from pathmover import (PathMover, MoveDetails, ReplicaExchange,
                       EnsembleHopMover)
from trajectory import Sample
from sample import SampleSet

class Calculation(object):

    calc_name = "Calculation"

    def __init__(self, storage, engine=None):
        self.storage = storage
        self.engine = engine

    def set_replicas(self, samples):
        self.globalstate = SampleSet(samples)

    def run(self, nsteps):
        print "Running an empty calculation? Try a subclass, maybe!"

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
        hopper = EnsembleHopMover(bias=self.bias,
                                  ensembles=[ensemble_from, ensemble_to],
                                  replicas=top_rep)

        shoot_samp = shooter.move(init_sample_set)
        init_sample_set.apply_samples(shoot_samp)
        hop_samp = hopper.move(init_sample_set)
        init_sample_set.apply_samples(hop_samp)

        # bring all the metadata from the submoves into our details
        details.__dict__.update(shoot_samp.details.__dict__)
        details.__dict__.update(hop_samp.details.__dict__)

        # set the rest of the details to their correct values
        details.inputs = details_inputs
        details.mover = details_mover
        details.replica = details_replica

        # the move will be accepted if the shooting move is accepted, no
        # matter what
        details.success = shoot_samp.details.success
        # the move is acceptable if the hopping is acceptable?
        details.accepted = hop_samp.details.accepted

        # result trajectory is whatever came out of hop_samp
        details.result = hop_samp.details.result

        setattr(details, 'start_replica', details.replica)
        if hop_samp.details.success == True:
            setattr(details, 'result_replica', details.replica+1)
        else:
            setattr(details, 'result_replica', details.replica)
        sample = Sample(replica=details.result_replica,
                        ensemble=details.result_ensemble,
                        trajectory=details.result,
                        details=details)

        #print "Success:", sample.details.success
        #print sample.trajectory, sample.details.result

        return sample


class Bootstrapping(Calculation):
    """The ensembles for the Bootstrapping calculation must be one ensemble
    set, in increasing order."""

    calc_name = "Bootstrapping"

    def __init__(self, storage, engine=None, movers=None, trajectory=None,
                 ensembles=None):
        super(Bootstrapping, self).__init__(storage, engine)
        self.ensembles = ensembles
        sample = Sample(replica=0, trajectory=trajectory, 
                        ensemble=self.ensembles[0], details=MoveDetails())
        self.globalstate = SampleSet([sample])
        if movers is None:
            pass # TODO: implement defaults: one per ensemble, uniform sel
        else:
            self.movers = movers

    def run(self, nsteps):
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
            print step_num, "Ensemble:", ens_num, "  failsteps=", failsteps
            old_rep = max(self.globalstate.replica_list())
            sample = bootstrapmove.move(self.globalstate)
            self.globalstate.apply_samples(sample, step=step_num)

            if sample.replica == old_rep:
                failsteps += 1
            else:
                failsteps = 0
                ens_num += 1

            if self.storage is not None:
                self.globalstate.save_samples(self.storage)
            step_num += 1

        for sample in self.globalstate:
            assert sample.ensemble(sample.trajectory) == True, "WTF?"

