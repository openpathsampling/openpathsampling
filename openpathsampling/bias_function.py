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

    def __add__(self, other):
        # the following craziness is to get the ensembles listed in the
        # right order: self-only, other-only, both; all preserving the order
        # in the original (preferring self's order for things in both)
        self_ensembles = self.ensembles_to_ids.keys()
        other_ensembles = other.ensembles_to_ids.keys()
        self_only = set(self_ensembles) - set(other_ensembles)
        other_only = set(other_ensembles) - set(self_ensembles)
        both = set(self_ensembles) & set(other_ensembles)
        all_ensembles = [e for e in self_ensembles if e in self_only]
        all_ensembles += [e for e in other_ensembles if e in other_only]
        all_ensembles += [e for e in self_ensembles if e in both]

        # set up the structures for initialization of the return
        ensembles_to_ids = {ens : all_ensembles.index(ens) 
                            for ens in all_ensembles}
        dataframe = pd.DataFrame(index=ensembles_to_ids.values(),
                                 columns=ensembles_to_ids.values())

        # fill the dataframe
        # to do this, we split the matrix into regions based on self_only,
        # other_only, and both. We have to check that we don't have two
        # versions of the 
        for ens_from in self_only:
            orig_from_id = self.ensembles_to_ids[ens_from]
            for ens_to in self_only | both:
                orig_to_id = self.ensembles_to_ids[ens_to]
                orig_value = self.dataframe.loc[orig_from_id, orig_to_id]
                dataframe.set_value(index=ensembles_to_ids[ens_from],
                                    col=ensembles_to_ids[ens_to],
                                    value=orig_value)

        for ens_from in other_only:
            orig_from_id = other.ensembles_to_ids[ens_from]
            for ens_to in other_only | both:
                orig_to_id = other.ensembles_to_ids[ens_to]
                orig_value = other.dataframe.loc[orig_from_id, orig_to_id]
                dataframe.set_value(index=ensembles_to_ids[ens_from],
                                    col=ensembles_to_ids[ens_to],
                                    value=orig_value)

        for ens_from in both:
            self_from_id = self.ensembles_to_ids[ens_from]
            other_from_id = other.ensembles_to_ids[ens_from]
            for ens_to in self_only:
                to_id = self.ensembles_to_ids[ens_to]
                orig_value = self.dataframe.loc[self_from_id, to_id]
                dataframe.set_value(index=ensembles_to_ids[ens_from],
                                    col=ensembles_to_ids[ens_to],
                                    value=orig_value)
            for ens_to in other_only:
                to_id = other.ensembles_to_ids[ens_to]
                orig_value = other.dataframe.loc[other_from_id, to_id]
                dataframe.set_value(index=ensembles_to_ids[ens_from],
                                    col=ensembles_to_ids[ens_to],
                                    value=orig_value)
            for ens_to in both:
                self_to_id = self.ensembles_to_ids[ens_to]
                other_to_id = self.ensembles_to_ids[ens_to]
                try:
                    self_value = self.dataframe.loc[self_from_id,
                                                    self_to_id]
                except KeyError:
                    # expected if this is only in other
                    self_value = None
                try:
                    other_value = other.dataframe.loc[other_from_id,
                                                      other_to_id]
                except KeyError:
                    # expected if this is only in self
                    other_value = None
                if self_value is not None and other_value is not None:
                    if self_value != other_value:
                        msg = "Biases have differing value for same entry:"
                        msg = " {:} != {:}".format(self_value, other_value)
                        raise ValueError(msg)
                value = self_value if self_value is not None else other_value
                if value is not None:
                    dataframe.set_value(index=ensembles_to_ids[ens_from],
                                        col=ensembles_to_ids[ens_to],
                                        value=value)
                # if both self_value and other_value are None, df unchanged

        return BiasEnsembleTable(dataframe, ensembles_to_ids)


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


