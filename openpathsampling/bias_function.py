from openpathsampling.netcdfplus import StorableNamedObject
import pandas as pd

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
    def __init__(self, dataframe, ensembles_to_ids):
        super(BiasEnsembleTable, self).__init__()
        self.dataframe = dataframe
        self.ensembles_to_ids = ensembles_to_ids

    @classmethod
    def ratios_from_dictionary(cls, ratio_dictionary):
        ensembles_to_ids = {e : ratio_dictionary.keys().index(e) 
                            for e in ratio_dictionary.keys()}
        id_based_df = pd.DataFrame(index=ensembles_to_ids.values(),
                                   columns=ensembles_to_ids.values())
        for e_from in ratio_dictionary:
            w_from = ratio_dictionary[e_from]
            id_from = ensembles_to_ids[e_from]
            for e_to in ratio_dictionary:
                id_to = ensembles_to_ids[e_to]
                w_to = ratio_dictionary[e_to]
                id_based_df.set_value(id_from, id_to, w_from/w_to)
        return BiasEnsembleTable(id_based_df, ensembles_to_ids)


    def bias_value(self, from_ensemble, to_ensemble):
        from_id = self.ensembles_to_ids[from_ensemble]
        to_id = self.ensembles_to_ids[to_ensemble]
        return self.dataframe.loc[from_id, to_id]

    @property
    def df(self):
        ids_to_ensembles = {self.ensembles_to_ids[e] : e
                            for e in self.ensembles_to_ids}
        df = self.dataframe.copy()
        df.index = [ids_to_ensembles[i].name for i in df.index]
        df.columns = [ids_to_ensembles[i].name for i in df.columns]
        return df.sort_index(axis=0).sort_index(axis=1)

    def bias_probability(self, from_ensemble, to_ensemble):
        # I have no idea why this needs ro be wrapped in a float to be
        # recognized as a float, but apparently it does
        return float(min(1.0, self.bias_value(from_ensemble, to_ensemble)))

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


