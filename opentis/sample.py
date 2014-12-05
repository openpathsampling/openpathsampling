
import random

class SampleSet(object):
    def __init__(self, samples):
        self.samples = samples
        self.ensemble_dict = {}
        self.replica_dict = {}
        for sample in samples:
            self.add(sample)

    def __getitem__(self, key):
        if isinstance(key, Ensemble):
            return self.ensemble_dict[ensemble]
        else:
            return self.replica_dict[replica]

    def __setitem__(self, key, value):
        try:
            dead_to_me = self[key]
        except KeyError:
            dead_to_me = None
        if dead_to_me is not None:
            del self[dead_to_me]
        self[key] = value

    def __delitem__(self, sample):
        del self.ensemble_dict[sample]
        del self.replica_dict[sample]
        del self.samples[sample]

    def __iter__(self):
        for sample in self.samples:
            yield sample

    def __len__(self):
        return len(self.samples)

    def __contains__(self, item):
        return (item in self.samples)

    def all_from_ensemble(ensemble):
        return self.ensemble_dict[ensemble]

    def all_from_replica(replica):
        return self.replica_dict[replica]

    def add(sample):
        try:
            self.ensemble_dict[sample.ensemble].append(sample)
        except KeyError:
            self.ensemble_dict[sample.ensemble] = [sample]
        try:
            self.replica_dict[sample.replica].append(sample)
        except KeyError:
            self.replica_dict[sample.replica] = [sample]

