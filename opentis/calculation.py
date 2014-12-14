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
        super(BootstrapPromotionMove, self).__init__(ensembles, replicas)
        self.shooters = shooters
        self.bias = bias

    def move(self, globalstate):
        # the tricky part here is that, if the hop is allowed, we only want
        # to report the sample in the new ensemble. The way we do this is by
        # treating each bootstrap move as a combination of 3 moves: it
        # always starts with a shooting move and a replica hop, and then, if
        # the hop was successful, a replica ID change move
        top_rep = max(globalstate.replica_list())
        ensemble_from = self.ensembles[top_rep]
        ensemble_to = self.ensembles[top_rep+1]
        old_sample = globalstate[top_rep]

        details = MoveDetails()
        details.inputs = [old_sample.trajectory]
        details.mover = self
        details.replica = top_rep

        init_sample_set = SampleSet([old_sample])

        shooter = self.shooters[top_rep]
        hopper = EnsembleHopMover(bias=self.bias,
                                  ensembles=[ensemble_from, ensemble_to],
                                  replicas=top_rep)
        #repchanger = ReplicaIDChange(new_replicas={top_rep : top_rep + 1}
                                     #old_samples={top_rep : old_sample},
                                     #ensembles=ensemble_from,
                                     #replcias=top_rep
                                    #)

        shoot_samp = shooter.move(init_sample_set)
        init_sample_set = init_sample_set.apply_samples(shoot_samp)
        hop_samp = hopper.move(init_sample_set)
        init_sample_set = init_sample_set.apply_samples(hop_samp)
        # TODO: this should be made into a function somewhere
        details.take_attributes_from(shoot_samp)
        details.take_attributes_from(hop_samp)

        details.success = shoot_samp.success
        details.accepted = hop_samp.accepted
        details.acceptance = shoot_samp.acceptance*hop_samp.acceptance
        details.result = hop_samp.result

        setattr(details, 'start_replica', sample.replica)
        if hop_samp.details.success == True:
            setattr(details, 'result_replica', sample.replica+1)
        else:
            setattr(details, 'result_replica', sample.replica)
        sample.replica = sample.details.result_replica
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
                        ensemble=self.ensembles[0])
        self.globalstate = SampleSet([sample])
        if movers is None:
            pass # TODO: implement defaults
        else:
            self.movers = movers

    def run(self, nsteps):

        bootstrapmove = BootstrapPromotionMove()

        ens_num = 0
        failsteps = 0
        # if we fail nsteps times in a row, kill the job

        while ens_num < len(self.ensembles) - 1 and failsteps < nsteps:
            # do a shooting move
            print self.globalstate.ensemble_dict
            print self.globalstate.replica_dict
            shoot_sample = self.movers[ens_num].move(self.globalstate)
            # bootstrap: if we crossed the next interface, promote
            hop_sample = bootstrapmove.move(SampleSet([tmp_sample]))

            if hop_sample.details.success == False:
                failsteps += 1
                to_apply = [shoot_sample, hop_sample]
            else:
                to_apply = [hop_sample]
            # apply the result to our current status
            self.globalstate = self.globalstate.apply_samples(to_apply)

            if self.storage is not None:
                self.globalstate.save_samples(self.storage)

        print "Done with new version"
        while ens_num < len(self.ensembles) - 1 and failsteps < nsteps:
#            print "Trying move in ensemble", ens_num
            # Generate Samples
            samples = [self.movers[ens_idx].move(self.globalstate[ens_idx]) 
                       for ens_idx in range(ens_num, ens_num + 1)]

            # Generate new globalstate using only the one sample
            globalstate = self.globalstate.apply_samples(samples)

            # Now save all samples
            if self.storage is not None:
                self.globalstate.save_samples(self.storage)

            # update to new globalstate
            self.globalstate = globalstate

            if ens_num < self.globalstate.size:
                # We can try to switch to the next ensemble

                # Check if the new trajectory (still in the old ensemble) would
                # fir in the next ensemble
                if self.globalstate.ensembles[ens_num + 1](self.globalstate[ens_num]):
                    # Yes, so apply the BootStrapMove and generate a new sample in the next ensemble
                    sample = bootstrapmove.move(self.globalstate[ens_num], 
                                                self.globalstate.ensembles[ens_num + 1])
                    globalstate = self.globalstate.apply_samples([sample])

                    # Now save all samples
                    self.globalstate.save_samples(self.storage)

                    # update to new globalstate
                    self.globalstate = globalstate

                    failsteps = 0
                    ens_num += 1

            failsteps += 1
