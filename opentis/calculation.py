
import os
import sys


class Calculation(object):

    calc_name = "Calculation"

    def __init__(self, storage, simulator=None, ensembles=None, 
                 movers=None):
        self.storage = storage
        self.simulator = simulator
        self.ensembles = ensembles
        self.movers = movers

    # TODO: this should be a property
    def set_replicas(replicas):
        self.replicas = replicas
    
    def run(self, nsteps):
        print "Running an empty calculation? Try a subclass, maybe!"


class Bootstrapping(Calculation):
    """The ensembles for the Bootstrapping calculation must be one ensemble
    set, in increasing order."""

    calc_name = "Bootstrapping"

    def run(self, nsteps):
        ens_num = 0
        failsteps = 0
        # if we fail nsteps times in a row, kill the job
        while ens_num < len(self.ensembles) and failsteps < nsteps: 
            print "Trying move in ensemble", ens_num
            sample = self.movers[ens_num].move(self.replicas[0])
            self.storage.sample.save(sample)
            self.replicas[0] = sample.trajectory
            # formally, I should probably create a move to promote the
            # trajectory to another ensemble and save that sample;
            # practically, I don't care
            add_to_ens = 0 if ens_num == len(self.ensembles) - 1 else 1
            if self.ensembles[ens_num+add_to_ens](sample.trajectory):
                failsteps = 0
                ens_num += 1
            failsteps += 1
