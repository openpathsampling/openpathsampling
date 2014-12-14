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
    def __init__(self, bias=None, ensembles=None, replicas='all'):
        super(BootstrapPromotionMove, self).__init__(ensembles, replicas)
        ens_neighbors = [ [a, b] for (a, b) in zip(ensembles, ensembles[1:]) ]
        self.hopper = EnsembleHopMover(bias, ens_neighbors, replicas)
        self.rep_id_dict = {}
        for ens in self.ensembles:
            self.rep_id_dict[ens] = self.ensembles.index(ens)

    def move(self, globalstate):
        intermed_samp = self.hopper.move(globalstate)

        details = MoveDetails()
        details.inputs = [trajectory]
        details.mover = self
        details.final = trajectory
        details.success = True
        details.acceptance = 1.0
        details.result = trajectory

        sample = Sample(
            trajectory=details.result,
            mover=self,
            ensemble=ensemble,
            details=details
        )

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

        swapmove = ReplicaExchange()
        bootstrapmove = BootstrapEnsembleChangeMove()

        ens_num = 0
        failsteps = 0
        # if we fail nsteps times in a row, kill the job
        while ens_num < len(self.globalstate) - 1 and failsteps < nsteps:
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
