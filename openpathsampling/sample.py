import random

import openpathsampling as paths

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

        JHP: I changed it to be immutable! This makes it slower but better
        tractable. Just a try

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

    # TODO: Can a sample be several times in an SampleSet? Why not?

    def __init__(self, samples = None, predecessor = None, accepted=True):
        self.predecessor = predecessor
        self.accepted = accepted
        self._samples = []
        self._final = None
        self._ensemble_dict = {}
        self._replica_dict = {}
        if samples is None:
            samples = []
        if type(samples) is Sample:
            samples = [samples]

        for sample in samples:
            if sample not in self._samples:
                self._append(sample)


    def _append(self, sample):
        self._samples.append(sample)
        try:
            self._ensemble_dict[sample.ensemble].append(sample)
        except KeyError:
            self._ensemble_dict[sample.ensemble] = [sample]
        try:
            self._replica_dict[sample.replica].append(sample)
        except KeyError:
            self._replica_dict[sample.replica] = [sample]

    @property
    def samples(self):
        return self._samples

    @property
    def final(self):
        if self._final is None:
            self._final = self + SampleSet(self._samples, self)

    @property
    def ensemble_dict(self):
        return self._ensemble_dict

    def replica_dict(self):
        return self._replica_dict

    def __plus__(self, other):
        if other.predecessor is self:
            newset = self.copy()
            for sample in other._samples:
                if sample not in self._samples:
                    self._append(sample)

            return newset
        else:
            raise ValueError('Incompatible SampleSets')

    def copy(self):
        new_set = SampleSet(self, accepted=self.accepted)

        return new_set

    def __getitem__(self, key):
        if isinstance(key, paths.Ensemble):
            return random.choice(self._ensemble_dict[key])
        else:
            return random.choice(self._replica_dict[key])

    def __setitem__(self, key, value):
        # In place substitution is not allowed for immutable items
        raise ValueError('no setitem for immutable lists!')

    def replace(this, key, value):
        # first, we check whether the key matches the sample: if no, KeyError
        if isinstance(key, paths.Ensemble):
            if key != value.ensemble:
                raise SampleKeyError(key, value, value.ensemble)
        else:
            if key != value.replica:
                raise SampleKeyError(key, value, value.replica)

        self = this.copy()

        if value in self.samples:
            # if value is already in this, we don't need to do anything
            return this
        # Setting works by replacing one with the same key. We pick one with
        # this key at random (using __getitem__), delete it, and then append
        # the new guy. If nothing exists with the desired key, this is the
        # same as append.
        try:
            dead_to_me = self[key]
        except KeyError:
            dead_to_me = None
        if dead_to_me is not None:
            self = self.delete(dead_to_me)

        self = self.append(value)
        return self

    def __delitem__(self, sample):
        # Same here. No delete for immutable lists
        # We could return a copy and redirect self to the copy, but
        # this could be confused since writing a[3] = b assumes that
        # a will not be replaced by a copy!
        # Also the underlying sample object can be changed!

        raise ValueError('no in place delete for immutable items')

    def delete(this, sample):
        self = this.copy()

        self._ensemble_dict[sample.ensemble].remove(sample)
        self._replica_dict[sample.replica].remove(sample)
        if len(self._ensemble_dict[sample.ensemble]) == 0:
            del self._ensemble_dict[sample.ensemble]
        if len(self._replica_dict[sample.replica]) == 0:
            del self._replica_dict[sample.replica]
        self._samples.remove(sample)

        return self

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
            return self._ensemble_dict[ensemble]
        except KeyError:
            return []

    def all_from_replica(self, replica):
        try:
            return self._replica_dict[replica]
        except KeyError:
            return []

    def append(self, sample):
        #  This works and will actually return a new object

        if sample in self._samples:
            # question: would it make sense to raise an error here? can't
            # have more than one copy of the same sample, but should we
            # ignore it silently or complain?
            return

        this = self.copy()

        this._samples.append(sample)
        try:
            this._ensemble_dict[sample.ensemble].append(sample)
        except KeyError:
            this._ensemble_dict[sample.ensemble] = [sample]
        try:
            this._replica_dict[sample.replica].append(sample)
        except KeyError:
            this._replica_dict[sample.replica] = [sample]

        return this

    def get_ensemble_dict(self):
        """
        Returns the dictionary of ensembles and their samples but not cached
        :return:
        """
        ensembles = set([sample.ensemble for sample in self._samples])
        return { sample.ensemble : [sample for sample in self._samples if sample.ensemble is ensemble] for ensemble in ensembles}

    def get_replica_dict(self):
        """
        Returns the dictionary of replica and their samples but not cached
        :return:
        """
        replicas = set([sample.replica for sample in self._samples])
        return { sample.ensemble : [sample for sample in self._samples if sample.replica is replica] for replica in replicas}


    def extend(this, samples):
        # note that this works whether the parameter samples is a list of
        # samples or a SampleSet!

        self = this.copy()
        try:
            for sample in samples:
                self = self.append(sample)
        except TypeError:
            # also acts as .append() if given a single sample
            self.append(samples)

        return self

    def apply(self, samples=None, accepted=None, move=None):
        '''Updates the SampleSet based on a list of samples, by setting them
        by replica in the order given in the argument list.'''
        if type(samples) is Sample:
            samples = [samples]

        if accepted is None:
            accepted = self.accepted

        if samples is None:
            samples = []

        newset = SampleSet(self, predecessor=self.predecessor, accepted=accepted)

        for sample in samples:
            # TODO: should time be a property of Sample or SampleSet?
            newset = newset.replace(sample.replica, sample)

        return newset

    def replica_list(self):
        '''Returns the list of replicas IDs in this SampleSet'''
        return self._replica_dict.keys()

    def ensemble_list(self):
        '''Returns the list of ensembles in this SampleSet'''
        return self._ensemble_dict.keys()
            
    def save_samples(self, storage):
        """
        Save all samples in the current GlobalState object. This should be
        called after a move has generated a new object since then all
        samples will get a timestamp that is associated with this

        This should not be necessary since storage of the SampleSet should
        trigger saving all samples

        Parameters
        ==========
        storage : Storage()
            the underlying netcdf file to be used for storage
        """
        map(storage.sample.save, self.samples)

    def consistency_check(self):
        '''This is mainly a sanity check for use in testing, but might be
        good to run (rarely) in the code until we're sure the tests cover
        all use cases.'''
        # check that we have the same number of samples in everything
        nsamps_ens = 0
        for ens in self._ensemble_dict.keys():
            nsamps_ens += len(self._ensemble_dict[ens])
        nsamps_rep = 0
        for rep in self._replica_dict.keys():
            nsamps_rep += len(self._replica_dict[rep])
        nsamps = len(self.samples)
        assert nsamps==nsamps_ens, \
                "nsamps != nsamps_ens : %d != %d" % (nsamps, nsamps_ens)
        assert nsamps==nsamps_rep, \
                "nsamps != nsamps_rep : %d != %d" % (nsamps, nsamps_rep)

        # if we have the same number of samples, then we check that each
        # sample in samples is in each of the dictionaries
        for samp in self.samples:
            assert samp in self._ensemble_dict[samp.ensemble], \
                    "Sample not in ensemble_dict! %r %r" % (samp, self._ensemble_dict)
            assert samp in self._replica_dict[samp.replica], \
                    "Sample not in replica_dict! %r %r" % (samp, self._replica_dict)

        # finally, check to be sure that thre are no duplicates in
        # self.samples; this completes the consistency check
        for samp in self.samples:
            assert self.samples.count(samp) == 1, \
                    "More than one instance of %r!" % samp

class EmptySampleSet(SampleSet):
    def __init__(self):

class InitialSampleSet(SampleSet):
    """
    A SampleSet with empty initial referenced sampleset
    """
    def __init__(self, samples=None, accepted=True):
        self.predecessor =


class SequentialSampleSet(SampleSet):
    def __init__(self, sets):


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

    def __init__(self, replica=None, trajectory=None, ensemble=None, details=None, step=-1):
        self.replica = replica
        self.ensemble = ensemble
        self.trajectory = trajectory
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

    @staticmethod
    def set_time(step, samples):
        for sample in samples:
            sample.step = step
