from openpathsampling.netcdfplus import StorableNamedObject

# NOTE: the biases that return here can still be more than 1. This is
# correct. You only fix them to min(1, value) in the Metropolis acceptance
# criterion.

class BiasFunction(StorableNamedObject):
    """Generic bias functions. Everything inherits from here. Abstract."""
    def probability_old_to_new(self, sampleset, change):
        raise NotImplementedError

    def probability_new_to_old(self, sampleset, change):
        raise NotImplementedError

    def get_new_old(self, sampleset, change):
        """Associates changed and original samples.

        Returns tuple (replica, new_sample, old_sample)
        """
        # TODO: this should move to sample.py (or pmc.py?) as a free function
        return [(s.replica, s, sampleset[s.replica]) for s in change.samples]


class BiasLookupFunction(BiasFunction):
    """Abstract class for lookup function based bias functions.
    """
    def __init__(self, bias_lookup, sample_reducer):
        super(BiasLookupFunction, self).__init__()
        self.bias_lookup = bias_lookup
        self.sample_reducer = sample_reducer


class BiasEnsembleTable(BiasFunction):
    # TODO: bias seems kind of fixed to Metropolis acceptance criterion --
    # is that okay elsewhere?
    def __init__(self, bias_table):
        super(BiasEnsembleTable, self).__init__()
        self.bias_table = bias_table

    def bias_value(self, from_ensemble, to_ensemble):
        from_w = self.bias_table[from_ensemble]
        to_w = self.bias_table[to_ensemble]
        return from_w / to_w

    def bias_probability(self, from_ensemble, to_ensemble):
        return min(1.0, self.bias_value(from_ensemble, to_ensemble))

    def probability_old_to_new(self, sampleset, change):
        new_old = self.get_new_old(sampleset, change)
        prob = 1.0
        for diff in new_old:
            new = diff[1]
            old = diff[2]
            prob *= self.bias_value(old.ensemble, new.ensemble)
        return min(1.0, prob)

    def probability_new_to_old(self, sampleset, change):
        new_old = self.get_new_old(sampleset, change)
        prob = 1.0
        for diff in new_old:
            new = diff[1]
            old = diff[2]
            prob *= self.bias_value(new.ensemble, old.ensemble)
        return min(1.0, prob)

