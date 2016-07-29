#!/usr/bin/env python

# wham.py
# The Weighted Histogram Analysis Method (WHAM) for combining histograms
# from several files.
import pandas as pd
import numpy as np

import logging
logger = logging.getLogger(__name__)

class WHAM(object):
    """
    Weighted Histogram Analysis Method

    Note
    ----
    Several parts of the docstrings mention F&S, which is intended to refer
    the reader to reference [1]_, in particular pages 184-187 in the 2nd
    edition (section called "Self-Consistent Histogram Method").

    Reference
    ---------
    .. [1] Daan Frenkel and Berend Smit. Understanding Molecular Simulation:
       From Algorithms to Applications. 2nd Edition. 2002.
    """

    def __init__(self, tol=1e-10, max_iter=1000000, cutoff=0.05):
        """
        Initialize (empty) WHAM calculation object.

        Parameters
        ----------
        tol : float
            tolerance for convergence
        max_iter : int
            maximum number of iterations
        cutoff : float
            windowing cutoff

        Attributes
        ----------
        sample_every : int
            frequency (in iterations) to report debug information

        """
        self.tol = tol
        self.max_iter = max_iter
        self.cutoff = cutoff

        self.sample_every = max_iter + 1
        self._float_format = "10.8"

    @property
    def float_format(self):
        """Float output format. Example: 10.8 (default)"""
        return lambda x : "{:" + self._float_format + "f}".format(x)

    @float_format.setter
    def float_format(self, value):
        self._float_format = value

    def load_files(self,fnames):  # pragma: no cover
        """Load a file or files into the internal structures.

        Requires either pandas or something else with pandas-like read_table
        and concat functions.

        Parameters
        ----------
        fnames : string or list of string
            file(s) to read in; each gives a column of the dataframe, with
            common indexes (keys)

        Returns
        -------
        pd.DataFrame
            dataframe with each input histogram in a column
        """
        # TODO: add some validation of formatting
        frames = []
        try:
            for fname in fnames:
                frames.append(pd.read_table(fname, index_col=0, sep=" ",
                                            usecols=[0,1], header=None))
        except TypeError:
            frames.append(pd.read_table(fnames, index_col=0, sep=" ",
                                        usecols=[0,1], header=None))
            fnames = [fnames]
        df = pd.concat(frames, axis=1)
        df.columns=fnames
        return df


    def prep_reverse_cumulative(self, df, cutoff=None, tol=None):
        """
        Created cleaned dataframe for further analysis.

        This version assumes that the input is the result of a reversed
        cumulative histogram. That means that it cleans leading input where
        the initial 
        

        Returns
        -------
        cleaned_df
        """
        if cutoff is None:
            cutoff = self.cutoff
        if tol is None:
            tol = self.tol

        # clear things that don't pass the cutoff
        hist_max = df.max(axis=0)
        raw_cutoff = cutoff*hist_max
        cleaned_df = df.apply(
            lambda s : [x if x > raw_cutoff[s.name] else 0.0 for x in s]
        )

        # clear duplicates of leading values
        test_f = lambda val1, val2, val_max : (
            abs(val1 - val2) > tol or abs(val1 - val_max) > tol
        )
        cleaned_df = cleaned_df.apply(
            lambda s : [
                s.iloc[i] if test_f(s.iloc[i], s.iloc[i+1], s.max()) else 0.0
                for i in range(len(s)-1)
            ] + [s.iloc[-1]]
        )
        return cleaned_df


    def unweighting_tis(self, cleaned_df):
        """
        Calculates the "unweighting" values for each histogram.

        In TIS, this is just 1 if there is a non-zero entry in cleaned_df,
        and 0 otherwise.
        """
        unweighting = cleaned_df.copy().applymap(
            lambda x : 1.0 if x > 0.0 else 0.0
        )
        return unweighting


    def sum_k_Hk_Q(self, cleaned_df):
        """
        Returns
        -------
        sum_k_Hk_Q 
            called `sum_hist` in other codes, or :math:`\sum_k H_k(Q)` in
            F&S. This is the sum over histograms of values for a given
            histogram bin.
        """
        return cleaned_df.sum(axis=1)


    def n_entries(self, cleaned_df):
        """
        n_entries : 
            the list of counts of entries. In other codes, this is `nt`. In
            F&S, this is :math:`M_k`.
        """
        return cleaned_df.sum(axis=0)


    def weighted_counts_tis(self, unweighting, n_entries):
        """
        Returns
        -------
        pd.Panel :
            weighted counts matrix, with 
        """
        weighted_counts = unweighting.apply(lambda s : [x * n_entries[s.name]
                                                        for x in s])
        return weighted_counts


    def generate_lnZ(self, lnZ, unweighting, weighted_counts,
                            sum_k_Hk_Q, tol=None):
        """
        Parameters
        ----------
        lnZ : list-like, one per histogram
            initial guess for ln(Z_i) for each histogram i
        unweighting : matrix-like, 
            the unweighting matrix for each histogram point. For TIS, this
            is 1 if the (cleaned) DF has an entry; 0 otherwise. In F&S, this
            is :math:`\exp(-\\beta W_i)`.
        sum_k_Hk_Q :
            see 
        """
        if tol is None:
            tol = self.tol
        diff = self.tol + 1  # always start above the tolerance
        iteration = 0
        hists = weighted_counts.columns
        bins = weighted_counts.index
        # TODO: probably faster if we make wc this a sparse matrix
        wc = weighted_counts.as_matrix()
        unw = unweighting.as_matrix()
        lnZ_old = pd.Series(data=lnZ, index=hists)
        Z_new = pd.Series(index=hists)
        sum_k_Hk_byQ = sum_k_Hk_Q.as_matrix()
        while diff > tol and iteration < self.max_iter:
            Z_old = np.exp(lnZ_old)
            reciprocal_Z_old = (1.0 / Z_old).as_matrix()
            for i in range(len(hists)):
                hist_i = hists[i]
                Z_new_i = 0.0

                w_i = unw[:,i]
                numerator_byQ = np.multiply(w_i, sum_k_Hk_byQ)
                sum_over_Z_byQ = wc.dot(reciprocal_Z_old)
                addends_k = np.divide(numerator_byQ, sum_over_Z_byQ)
                Z_new_i = np.nansum(addends_k)

                # for val_q in range(len(bins)):
                    # local_w_i = unw[val_q, i] #unweighting.iloc[val_q, i]
                    # if local_w_i > 0:
                        # if statement allows us to skip if the weight is 0
                        # this will be \sum_j x_j M_j / Z_j^{(old)}
                        # sum_w_over_Z = np.dot(wc[val_q], reciprocal_Z_old)
                        # sum_w_over_Z = sum([
                            # (wc[val_q, hist_j] / Z_old.iloc[hist_j])
                            # for hist_j in range(len(hists))
                        # ])

                        # Z_new_i += (local_w_i * sum_k_Hk_Q.iloc[val_q] / sum_w_over_Z)
                Z_new[hist_i] = Z_new_i

            lnZ_new = np.log(Z_new)

            iteration += 1
            diff = self.get_diff(lnZ_old, lnZ_new, iteration)
            lnZ_old = lnZ_new - lnZ_new.iloc[0]

        logger.info("iterations=" + str(iteration) + " diff=" + str(diff))
        logger.info("       lnZ=" + str(lnZ_old))
        self.convergence = (iteration, diff)
        return lnZ_old


    def get_diff(self, lnZ_old, lnZ_new, iteration):
        # get error
        diff=0
        diff = sum(abs(lnZ_old - lnZ_new))
        # check status (mainly for debugging)
        if (iteration % self.sample_every == 0):  # pragma: no cover
            logger.debug("niteration = " + str(iteration))
            logger.debug("  diff = " + str(diff))
            logger.debug("   lnZ = " + str(lnZ_old))
            logger.debug("lnZnew = " + str(lnZ_new))
        return diff


    def output_histogram(self, lnZ, sum_k_Hk_Q, weighted_counts):
        Z = np.exp(lnZ)
        Z0_over_Zi = Z.iloc[0] / Z
        output = pd.Series(index=sum_k_Hk_Q.index, name="WHAM")
        for val in sum_k_Hk_Q.index:
            sum_w_over_Z = sum([
                weighted_counts.loc[val, hist_i] * Z0_over_Zi[hist_i]
                for hist_i in Z.index
            ])
            output[val] = sum_k_Hk_Q[val] / sum_w_over_Z

        return output

    @staticmethod
    def normalize_cumulative(series):
        return series/series.max()

    def guess_lnZ_crossing_probability(self, cleaned_df):
        df = cleaned_df.apply(lambda s : s/s.max())
        # pandas magic, see http://stackoverflow.com/questions/18327624
        first_nonzero = df.apply(lambda s: s[s == 1.0].index[0])
        df = df.loc[first_nonzero]
        guess_nextZ_over_Z = df.apply(lambda s: s.nlargest(2).iloc[1])
        guess_Z_data = [1.0]
        for i in range(len(guess_nextZ_over_Z)-1):
            guess_Z_data.append(guess_Z_data[-1] * guess_nextZ_over_Z.iloc[i])
        guess_lnZ = pd.Series(data=np.log(np.array(guess_Z_data)),
                              index=guess_nextZ_over_Z.index)
        return guess_lnZ

    def wham_bam_histogram(self, input_df):
        cleaned = self.prep_reverse_cumulative(input_df)
        guess = self.guess_lnZ_crossing_probability(cleaned)
        sum_k_Hk_Q = self.sum_k_Hk_Q(cleaned)
        n_entries = self.n_entries(cleaned)
        unweighting = self.unweighting_tis(cleaned)
        weighted_counts = self.weighted_counts_tis(unweighting,
                                                          n_entries)
        try:
            lnZ = self.generate_lnZ(guess, unweighting, weighted_counts,
                                    sum_k_Hk_Q)
        except IndexError as e:
            failmsg = "Does your input to WHAM have enough data?"
            if not e.args:
                e.args = [failmsg]
            else:
                arg0 = e.args[0]+"\n"+failmsg
                e.args = tuple([arg0] + list(e.args[1:]))
                raise e

        hist = self.output_histogram(lnZ, sum_k_Hk_Q, weighted_counts)
        return self.normalize_cumulative(hist)


def parsing(parseargs):  # pragma: no cover
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--tol", type="float", default=1e-12)
    parser.add_option("--max_iter", type="int", default=1000000)
    parser.add_option("--cutoff", type="float", default=0.05)
    parser.add_option("--pstats", type="string", default=None)
    parser.add_option("--float_format", type="string", default="10.8")
    opts, args = parser.parse_args(parseargs)
    return opts, args


import sys, os
if __name__ == "__main__":  # pragma: no cover
    opts, args = parsing(sys.argv[1:])
    wham = WHAM(tol=opts.tol, max_iter=opts.max_iter, cutoff=opts.cutoff)
    wham.float_format = opts.float_format
    df = wham.load_files(args)
    # wham.sample_every = 10

    if opts.pstats is not None:
        import cProfile
        import time
        start = time.time()
        cProfile.run("print wham.wham_bam_histogram(df)", opts.pstats)
        print time.time() - start
    else:
        wham_hist = wham.wham_bam_histogram(df)
        print wham_hist.to_string(header=False, 
                                  float_format=lambda x : "{:10.8f}".format(x))



