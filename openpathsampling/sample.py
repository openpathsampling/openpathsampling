import random

import openpathsampling as paths
import copy

class SampleKeyError(Exception):
    def __init__(self, key, sample, sample_key):
        self.key = key
        self.sample = sample
        self.sample_key = sample_key
        self.msg = (str(self.key) + " does not match " + str(self.sample_key)
                    + " from " + str(self.sample))

class SampleSet(object):
    '''
    SampleSet is essentially a list of samples, with a few conveniences.  It
    can be treated as a list of samples (using, e.g., .append), or as a
    dictionary of ensembles mapping to a list of samples, or as a dictionary
    of replica IDs to samples. Any type is allowed as a replica ID except
    Sample or Ensemble.

    The dictionaries ensemble_dict and replica_dict are conveniences which
    should be kept consistent by any method which modifies the container.
    They do not need to be stored.

    Note
    ----
        Current implementation is as an unordered set. Therefore we don't
        have some of the convenient tools in Python sequences (e.g.,
        slices). On the other hand, I'm not sure whether that is meaningful
        here.

    Attributes
    ----------
    samples : list of Sample
        The samples included in this set.
    ensemble_dict : dict
        A dictionary with Ensemble objects as keys and lists of Samples as
        values.
    replica_dict : dict
        A dictionary with replica IDs as keys and lists of Samples as values
    '''

    def __init__(self, samples, movepath=None):
        self.samples = []
        self.ensemble_dict = {}
        self.replica_dict = {}
        self.extend(samples)
        if movepath is None:
            self.movepath = paths.EmptyMovePath()
        else:
            self.movepath = movepath


    def __getitem__(self, key):
        if isinstance(key, paths.Ensemble):
            return random.choice(self.ensemble_dict[key])
        else:
            return random.choice(self.replica_dict[key])

    def __setitem__(self, key, value):
        # first, we check whether the key matches the sample: if no, KeyError
        if isinstance(key, paths.Ensemble):
            if key != value.ensemble:
                raise SampleKeyError(key, value, value.ensemble)
        else:
            if key != value.replica:
                raise SampleKeyError(key, value, value.replica)

        if value in self.samples:
            # if value is already in this, we don't need to do anything
            return
        # Setting works by replacing one with the same key. We pick one with
        # this key at random (using __getitem__), delete it, and then append
        # the new guy. If nothing exists with the desired key, this is the
        # same as append.
        try:
            dead_to_me = self[key]
        except KeyError:
            dead_to_me = None
        if dead_to_me is not None:
            del self[dead_to_me]
        self.append(value)

    def __delitem__(self, sample):
        self.ensemble_dict[sample.ensemble].remove(sample)
        self.replica_dict[sample.replica].remove(sample)
        if len(self.ensemble_dict[sample.ensemble]) == 0:
            del self.ensemble_dict[sample.ensemble]
        if len(self.replica_dict[sample.replica]) == 0:
            del self.replica_dict[sample.replica]
        self.samples.remove(sample)

    # TODO: add support for remove and pop

    def __iter__(self):
        for sample in self.samples:
            yield sample

    def __len__(self):
        return len(self.samples)

    def __contains__(self, item):
        return (item in self.samples)

    def all_from_ensemble(self, ensemble):
        try:
            return self.ensemble_dict[ensemble]
        except KeyError:
            return []

    def all_from_replica(self, replica):
        try:
            return self.replica_dict[replica]
        except KeyError:
            return []

    def append(self, sample):
        if sample in self.samples:
            # question: would it make sense to raise an error here? can't
            # have more than one copy of the same sample, but should we
            # ignore it silently or complain?
            return

        self.samples.append(sample)
        try:
            self.ensemble_dict[sample.ensemble].append(sample)
        except KeyError:
            self.ensemble_dict[sample.ensemble] = [sample]
        try:
            self.replica_dict[sample.replica].append(sample)
        except KeyError:
            self.replica_dict[sample.replica] = [sample]

    def extend(self, samples):
        # note that this works whether the parameter samples is a list of
        # samples or a SampleSet!
        try:
            for sample in samples:
                self.append(sample)
        except TypeError:
            # also acts as .append() if given a single sample
            self.append(samples)

    def apply_samples(self, samples, step=None, copy=True):
        '''Updates the SampleSet based on a list of samples, by setting them
        by replica in the order given in the argument list.'''
        if type(samples) is Sample:
            samples = [samples]
        if copy==True:
            newset = SampleSet(self)
        else:
            newset = self
        for sample in samples:
            # TODO: should time be a property of Sample or SampleSet?
            sample.step = step
            if sample.intermediate == False:
                newset[sample.replica] = sample
        return newset

    def apply_intermediates(self, samples, step=None, copy=True):
        '''Return updated SampleSet, including all intermediates.

        Useful in SequentialMovers.
        '''
        if type(samples) is Sample:
            samples = [samples]
        if copy==True:
            newset = SampleSet(self)
        else:
            newset = self
        for sample in samples:
            # TODO: should time be a property of Sample or SampleSet?
            sample.step = step
            newset[sample.replica] = sample
        return newset

    def replica_list(self):
        '''Returns the list of replicas IDs in this SampleSet'''
        return self.replica_dict.keys()

    def ensemble_list(self):
        '''Returns the list of ensembles in this SampleSet'''
        return self.ensemble_dict.keys()

    def save_samples(self, storage):
        """
        Save all samples in the current GlobalState object. This should be
        called after a move has generated a new object since then all
        samples will get a timestamp that is associated with this

        Parameters
        ==========
        storage : Storage()
            the underlying netcdf file to be used for storage
        """
        map(storage.sample.save, self.samples)

    def sanity_check(self):
        '''Checks that the sample trajectories satisfy their ensembles
        '''
        for sample in self:
            assert(sample.ensemble(sample.trajectory))

    def consistency_check(self):
        '''This is mainly a sanity check for use in testing, but might be
        good to run (rarely) in the code until we're sure the tests cover
        all use cases.'''
        # check that we have the same number of samples in everything
        nsamps_ens = 0
        for ens in self.ensemble_dict.keys():
            nsamps_ens += len(self.ensemble_dict[ens])
        nsamps_rep = 0
        for rep in self.replica_dict.keys():
            nsamps_rep += len(self.replica_dict[rep])
        nsamps = len(self.samples)
        assert nsamps==nsamps_ens, \
                "nsamps != nsamps_ens : %d != %d" % (nsamps, nsamps_ens)
        assert nsamps==nsamps_rep, \
                "nsamps != nsamps_rep : %d != %d" % (nsamps, nsamps_rep)

        # if we have the same number of samples, then we check that each
        # sample in samples is in each of the dictionaries
        for samp in self.samples:
            assert samp in self.ensemble_dict[samp.ensemble], \
                    "Sample not in ensemble_dict! %r %r" % (samp, self.ensemble_dict)
            assert samp in self.replica_dict[samp.replica], \
                    "Sample not in replica_dict! %r %r" % (samp, self.replica_dict)

        # finally, check to be sure that thre are no duplicates in
        # self.samples; this completes the consistency check
        for samp in self.samples:
            assert self.samples.count(samp) == 1, \
                    "More than one instance of %r!" % samp

    def __add__(self, other):
        """
        Add the move path to the Sample and return the new sample set
        """
        new_set = other.apply_to(self)
        new_set.movepath = other
        return new_set

    # @property
    # def ensemble_dict(self):
    #     if self._ensemble_dict is None:
    #         self._ensemble_dict = self._get_ensemble_dict()
    #
    #     return self._ensemble_dict
    #
    # def _get_ensemble_dict(self):
    #     """
    #     Returns the dictionary of ensembles and their samples but not cached
    #     :return:
    #     """
    #     ensembles = set([sample.ensemble for sample in self.samples])
    #     print ensembles
    #     return { sample.ensemble : [sample for sample in self.samples if sample.ensemble is ensemble] for ensemble in ensembles}
    #
    #
    # @property
    # def replica_dict(self):
    #     if self._replica_dict is None:
    #         self._replica_dict = self._get_replica_dict()
    #
    #     return self._replica_dict
    #
    # def _get_replica_dict(self):
    #     """
    #     Returns the dictionary of replica and their samples but not cached
    #     :return:
    #     """
    #     replicas = set([sample.replica for sample in self.samples])
    #     return { sample.replica : [sample for sample in self.samples if sample.replica is replica] for replica in replicas}
    #
    # def __plus__(self, other):
    #     if other.predecessor is self:
    #         newset = self.copy()
    #         for sample in other._samples:
    #             if sample not in self._samples:
    #                 self._append(sample)
    #
    #         return newset
    #     else:
    #         raise ValueError('Incompatible MovePaths')


class Sample(object):
    """
    A Sample represents a given "draw" from its ensemble, and is the return
    object from a PathMover. It and contains all information about the move,
    initial trajectories, new trajectories (both as references). 
    
    Since each Sample is a single representative of a single ensemble, each
    Sample consists of one replica ID, one trajectory, and one ensemble.
    This means that movers which generate more than one "draw" (often from
    different ensembles, e.g. replica exchange) will generate more than one
    Sample object.

    Attributes
    ----------
    replica : integer
        The replica ID to which this Sample applies
    trajectory : Trajectory
        The trajectory (path) for this sample
    ensemble : Ensemble
        The Ensemble this sample is drawn from
    details : MoveDetails
        Object 
    step : integer
        the Monte Carlo step number associated with this Sample
    """

    def __init__(self, replica=None, trajectory=None, ensemble=None, intermediate=False, details=None, step=-1):
        self.replica = replica
        self.ensemble = ensemble
        self.trajectory = trajectory
        self.intermediate = intermediate
        self.details = details
        self.step = step

    def __call__(self):
        return self.trajectory

    def __str__(self):
        mystr = "Step: "+str(self.step)+"\n"
        mystr += "Replica: "+str(self.replica)+"\n"
        mystr += "Trajectory: "+str(self.trajectory)+"\n"
        mystr += "Ensemble: "+repr(self.ensemble)+"\n"
        mystr += "Details: "+str(self.details)+"\n"
        return mystr

    def __repr__(self):
        return '<Sample @ ' + str(hex(id(self))) + '>'

    @staticmethod
    def set_time(step, samples):
        for sample in samples:
            sample.step = step

    def copy_reset(self):
        '''
        Copy of Sample with initialization move details.
        '''
        result = Sample(
            replica=self.replica,
            trajectory=self.trajectory,
            ensemble=self.ensemble,
            details=paths.MoveDetails.initialization(self)
        )
        return result


