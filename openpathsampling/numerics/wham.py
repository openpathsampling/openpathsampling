#!/usr/bin/env python

from __future__ import print_function

# wham.py
# The Weighted Histogram Analysis Method (WHAM) for combining histograms
# from several files.
import pandas as pd
import numpy as np
import sys
import logging
logger = logging.getLogger(__name__)


class WHAM(object):
    """
    Weighted Histogram Analysis Method

    Notes
    -----
        Several parts of the docstrings mention F&S, which is intended to
        refer the reader to reference [1]_, in particular pages 184-187 in
        the 2nd edition (section called "Self-Consistent Histogram Method").

        Other terminology: n_hists refers to the number of histograms,
        n_bins refers to the number of bins per histogram. Thus the input is
        a matrix of n_bins rows and n_hists columns.


    References
    ----------
    .. [1] Daan Frenkel and Berend Smit. Understanding Molecular Simulation:
       From Algorithms to Applications. 2nd Edition. 2002.

    Parameters
    ----------
    tol : float
        tolerance for convergence or equality. default 10e-10
    max_iter : int
        maximum number of iterations. Default 1000000
    cutoff : float
        windowing cutoff, as fraction of maximum value. Default 0.05

    Attributes
    ----------
    sample_every : int
        frequency (in iterations) to report debug information
    """
    def __init__(self, tol=1e-10, max_iter=1000000, cutoff=0.05,
                 interfaces=None):
        self.tol = tol
        self.max_iter = max_iter
        self.cutoff = cutoff
        self.interfaces = interfaces

        self.sample_every = max_iter + 1
        self._float_format = "10.8"
        self.lnZ = None

    @property
    def float_format(self):
        """Float output format. Example: 10.8 (default)"""
        return lambda x: "{:" + self._float_format + "f}".format(x)

    @float_format.setter
    def float_format(self, value):
        self._float_format = value

    def load_files(self, fnames):  # pragma: no cover
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
                                            usecols=[0, 1], header=None))
        except TypeError:
            frames.append(pd.read_table(fnames, index_col=0, sep=" ",
                                        usecols=[0, 1], header=None))
            fnames = [fnames]
        df = pd.concat(frames, axis=1)
        df.columns = fnames
        return df

    def prep_reverse_cumulative(self, df, cutoff=None, tol=None):
        """
        Created cleaned dataframe for further analysis.

        This version assumes that the input is the result of a reversed
        cumulative histogram. That means that it cleans leading input where
        the initial

        Parameters
        ----------
        df : pandas.DataFrame
            input dataframe
        cutoff : float
            windowing cutoff, as fraction of maximum value
        tol : float
            tolerance for two values being "equal"


        Returns
        -------
        pandas.DataFrame
            version of the dataframe with leading constant values removed
            and without anything under the cutoff. This is often referred to
            as the "cleaned input dataframe" in other functions
        """
        if cutoff is None:
            cutoff = self.cutoff
        if tol is None:
            tol = self.tol

        # clear things that don't pass the cutoff
        hist_max = df.max(axis=0)
        raw_cutoff = cutoff*hist_max
        cleaned_df = df.apply(
            lambda s: [x if x > raw_cutoff[s.name] else 0.0 for x in s]
        )

        if self.interfaces is not None:
            # use the interfaces values to set anything before that value to
            # zero
            if type(self.interfaces) is not pd.Series:
                self.interfaces = pd.Series(data=self.interfaces,
                                            index=df.columns)
            greater_almost_equal = lambda a, b: (a >= b or
                                                 abs(a - b) < 10e-10)
            cleaned_df = cleaned_df.apply(
                lambda s: [
                    (
                        s.iloc[i]
                        if greater_almost_equal(s.index[i],
                                                self.interfaces[s.name])
                        else 0.0
                    )
                    for i in range(len(s))
                ]
            )
        else:
            # clear duplicates of leading values
            test_f = lambda val1, val2, val_max: (
                abs(val1 - val2) > tol or abs(val1 - val_max) > tol
            )
            cleaned_df = cleaned_df.apply(
                lambda s: [
                    s.iloc[i]
                    if test_f(s.iloc[i], s.iloc[i+1], s.max()) else 0.0
                    for i in range(len(s)-1)
                ] + [s.iloc[-1]]
            )
        return cleaned_df

    def unweighting_tis(self, cleaned_df):
        """
        Calculates the "unweighting" values for each histogram.

        In TIS, this is just 1 if there is a non-zero entry in cleaned_df,
        and 0 otherwise.

        (NB: it isn't quite clear to me why this is a matrix instead of a
        vector, but the previous code accounted for a dependence on the
        bin of the histogram)

        Parameters
        ----------
        cleaned_df : pandas.DataFrame
            cleaned input dataframe

        Returns
        -------
        pandas.DataFrame
            unweighting values for the input dataframe
        """
        unweighting = cleaned_df.copy().applymap(
            lambda x: 1.0 if x > 0.0 else 0.0
        )
        return unweighting

    def sum_k_Hk_Q(self, cleaned_df):
        r"""Sum over histograms for each histogram bin. Length is n_bins

        Called ``sum_hist`` in other codes, or :math:`\sum_k H_k(Q)` in F&S.
        This is the sum over histograms of values for a given histogram bin.

        Parameters
        ----------
        cleaned_df : pandas.DataFrame
            cleaned input dataframe

        Returns
        -------
        pandas.Series
            sum over histograms for each bin (length is n_bins)
        """
        return cleaned_df.sum(axis=1)

    def n_entries(self, cleaned_df):
        """List of counts of entries per histogram. Length is n_hists

        The list of counts of entries. In other codes, this is `nt`. In F&S,
        this is :math:`M_k`.

        Parameters
        ----------
        cleaned_df : pandas.DataFrame
            cleaned input dataframe

        Returns
        -------
        pandas.Series
            List of counts of entries per histogram (length n_hists)
        """
        return cleaned_df.sum(axis=0)

    def weighted_counts_tis(self, unweighting, n_entries):
        """
        Product of unweighting and n_entries.

        In F&S, this is :math:`e^{-\\beta W_k} M_k`. This value appears as
        part of a sum in the WHAM iteration equation (F&S 2nd edition
        Eq. 7.3.10). The product can be pre-calculated, which is what we do
        here.

        As for this being a matrix (not a vector), see note on
        :meth:`.unweighting`.

        Returns
        -------
        pandas.DataFrame
            weighted counts matrix, size n_hists by n_dims
        """
        weighted_counts = unweighting.apply(lambda s: [x * n_entries[s.name]
                                                       for x in s])
        return weighted_counts

    def generate_lnZ(self, lnZ, unweighting, weighted_counts, sum_k_Hk_Q,
                     tol=None):
        r"""
        Perform the WHAM iteration to estimate ln(Z_i) for each histogram.

        Parameters
        ----------
        lnZ : pandas.Series, one per histogram (length n_hists)
            initial guess for ln(Z_i) for each histogram i
        unweighting : pandas.DataFrame, n_bins by n_hists
            the unweighting matrix for each histogram point. For TIS, this
            is 1 if the (cleaned) DF has an entry; 0 otherwise. In F&S, this
            is :math:`\exp(-\\beta W_i)`. See :meth:`.unweighting`.
        weighted_counts : pandas.DataFrame, n_bins by n_hists
            the weighted matrix multiplied by the counts per histogram. See
            :meth:`.weighted_counts_tis`.
        sum_k_Hk_Q : pandas.Series, one per bin (length n_bins)
            Sum over histograms for each histogram bin. See
            :meth:`.sum_k_Hk_Q`.
        tol : float
            convergence tolerance

        Returns
        -------
        pandas.Series
            the resulting WHAM calculation for ln(Z_i) for each histogram i
        """
        if tol is None:
            tol = self.tol
        diff = self.tol + 1  # always start above the tolerance
        iteration = 0
        hists = weighted_counts.columns
        # TODO: probably faster if we make wc this a sparse matrix
        wc = weighted_counts.values
        unw = unweighting.values
        lnZ_old = pd.Series(data=lnZ, index=hists)
        Z_new = pd.Series(index=hists, dtype='float64')
        sum_k_Hk_byQ = sum_k_Hk_Q.values
        while diff > tol and iteration < self.max_iter:
            Z_old = np.exp(lnZ_old)
            reciprocal_Z_old = (1.0 / Z_old).values
            for i in range(len(hists)):
                #############################################################
                # this is equation 7.3.10 in F&S
                # Z_i^{(new)} =
                #    \int \dd{Q} w_{i,Q}
                #    \times \frac{\sum_{j=1}^n H_j(Q)}
                #                {\sum_{k=1}^n w_{k,Q} M_k / Z_k^{(old)}}
                # where F&S explicitly use w_{i,Q} = e^{-\beta W_i}
                #
                # Matching terms from F&S to our variables:
                #   w_i = w_{i,Q} = $e^{-\beta W_i}$
                #       * this is a column (len == n_bins) from a matrix
                #         size n_hists \times n_bins
                #       * from "unweighting", which is Boltzmann in umbrella
                #         sampling (F&S), but 1 or 0 in TIS
                #       * see unweighting_tis: not entirely clear why this
                #         is a matrix, not an vector length n_hists
                #   sum_k_Hk_byQ = $\sum_{j=1}^n H_j(Q)$
                #       * this is a function of Q, thus len == n_bins
                #   wc = w_{k,Q} * M_k = $e^{-\beta W_k} M_k$
                #       * note that this is element-wise multiplication
                #       * matrix, size n_hists \times n_bins
                #   reciprocal_Z_old = $1/Z_k^{(old)}$
                #       * vector, len == n_hists
                #
                # Note that we do all of this with matrix/vector
                # multiplication, which numpy can do very quickly.
                #############################################################

                w_i = unw[:, i]

                # numerator: w_{i,Q} * sum_k_Hk_byQ
                numerator_byQ = np.multiply(w_i, sum_k_Hk_byQ)

                # denominator: wc * Z^{-1}
                sum_over_Z_byQ = wc.dot(reciprocal_Z_old)

                # divide each entry, and add them (integrate over Q in F&S)
                # we intentially allow invalid (0/0) to give NaN; gets
                # removed by using the np.nansum)
                with np.errstate(divide='ignore', invalid='ignore'):
                    addends_k = np.divide(numerator_byQ, sum_over_Z_byQ)
                Z_new[hists[i]] = np.nansum(addends_k)

            lnZ_new = np.log(Z_new)

            iteration += 1
            diff = self.get_diff(lnZ_old, lnZ_new, iteration)
            lnZ_old = lnZ_new - lnZ_new.iloc[0]

        logger.info("iterations=" + str(iteration) + " diff=" + str(diff))
        logger.info("       lnZ=" + str(lnZ_old))
        self.convergence = (iteration, diff)
        return lnZ_old

    def get_diff(self, lnZ_old, lnZ_new, iteration):
        """Calculate the difference for this iteration.

        Also outputs debug information, if desired.

        Parameters
        ----------
        lnZ_old : pandas.Series
            previous value of ln(Z_i)
        lnZ_new : pandas.Series
            new value of ln(Z_i)
        iteration : int
            iteration number

        Returns
        -------
        float
            difference between old and new to use for convergence testing
        """
        # get error
        diff = 0
        diff = sum(abs(lnZ_old - lnZ_new))
        # check status (mainly for debugging)
        if (iteration % self.sample_every == 0):  # pragma: no cover
            logger.debug("niteration = " + str(iteration))
            logger.debug("  diff = " + str(diff))
            logger.debug("   lnZ = " + str(lnZ_old))
            logger.debug("lnZnew = " + str(lnZ_new))
        return diff

    def output_histogram(self, lnZ, sum_k_Hk_Q, weighted_counts):
        """Recombine the data into a joined histogram

        Parameters
        ----------
        lnZ : pandas.Series, one per histogram (length n_hists)
            value of ln(Z_i) for each histogram i
        sum_k_Hk_Q : pandas.Series, one per bin (length n_bins)
            Sum over histograms for each histogram bin. See
            :meth:`.sum_k_Hk_Q`.
        weighted_counts : pandas.DataFrame, n_bins by n_hists
            the weighted matrix multiplied by the counts per histogram. See
            :meth:`.weighted_counts_tis`.

        Returns
        -------
        pandas.Series
            the WHAM-reweighted combined histogram, unnormalized
        """
        Z = np.exp(lnZ)
        Z0_over_Zi = Z.iloc[0] / Z
        output = pd.Series(index=sum_k_Hk_Q.index, name="WHAM",
                           dtype='float64')
        for val in sum_k_Hk_Q.index:
            sum_w_over_Z = sum([
                weighted_counts.loc[val, hist_i] * Z0_over_Zi[hist_i]
                for hist_i in Z.index
            ])
            # explicitly allow NaN results for simplcity (should only occur
            # when numerator and denominator are 0) ... this will leave NaNs
            # in the histogram in those locations; if all values of the
            # total histogram are NaN, that gets caught in the main
            # wham_bam_histogram routine
            with np.errstate(divide='ignore', invalid='ignore'):
                output[val] = sum_k_Hk_Q[val] / sum_w_over_Z

        return output

    @staticmethod
    def normalize_cumulative(series):
        """Normalize to maximum value"""
        return series/series.max()

    def guess_lnZ_crossing_probability(self, cleaned_df):
        """Guess ln(Z_i) based on crossing probabilities

        Parameters
        ----------
        cleaned_df : pandas.DataFrame
            cleaned input dataframe representing crossing probabilities

        Returns
        -------
        pandas.Series
            initial guess for values of ln(Z_i) for each histogram
        """
        df = cleaned_df.apply(lambda s: s/s.max())
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

    def check_cleaned_overlaps(self, cleaned_df):
        """
        Check that all the histograms have sufficient overlaps.

        Parameters
        ----------
        cleaned_df : pandas.DataFrame
            cleaned input dataframe representing crossing probabilities

        Raises
        ------
        RuntimeError
            if the input doesn't have enough of an overlap
        """
        for col_idx in range(1, len(cleaned_df.columns)):
            col = cleaned_df.columns[col_idx]
            prev_col = cleaned_df.columns[col_idx - 1]

            col_data = cleaned_df[col]
            prev_data = cleaned_df[prev_col]

            first_nonzero = col_data[col_data != 0.0].index[0]
            if not prev_data[first_nonzero] > 0.0:
                # use not and > to account for NaNs
                raise RuntimeError(
                    "Insufficient overlap to combine histograms.\n"
                    + "This is either due to poor sampling or bad "
                    + "interface placement.\n" + "Row: "
                    + str(first_nonzero) + "   Column: " + str(col) + "\n"
                    + str(cleaned_df))

    def wham_bam_histogram(self, input_df):
        """
        Perform the entire wham process.

        Parameters
        ----------
        input_df : pandas.DataFrame
            input dataframe

        Returns
        -------
        pandas.Series
            the WHAM-reweighted combined histogram, normalized so max value
            is 1
        """
        cleaned = self.prep_reverse_cumulative(input_df)
        self.check_cleaned_overlaps(cleaned)
        guess = self.guess_lnZ_crossing_probability(cleaned)
        sum_k_Hk_Q = self.sum_k_Hk_Q(cleaned)
        n_entries = self.n_entries(cleaned)
        unweighting = self.unweighting_tis(cleaned)
        weighted_counts = self.weighted_counts_tis(unweighting,
                                                   n_entries)
        try:
            lnZ = self.generate_lnZ(guess, unweighting, weighted_counts,
                                    sum_k_Hk_Q)
        except IndexError as e:  # pragma: no cover
            # I don't think this can happen any more, but leave it in case
            # (Now the check_cleaned_overlaps should catch this problem.)
            failmsg = "Does your input to WHAM have enough data?"
            if not e.args:
                e.args = [failmsg]
            else:
                arg0 = e.args[0]+"\n"+failmsg
                e.args = tuple([arg0] + list(e.args[1:]))
                raise e

        hist = self.output_histogram(lnZ, sum_k_Hk_Q, weighted_counts)
        result = self.normalize_cumulative(hist)
        self.lnZ = lnZ
        if sum(pd.isnull(result)) == len(result):  # pragma: no cover
            # last safety check
            raise RuntimeError("WHAM result is all NaN. Reason unknown.")
        return result


def parsing(parseargs):  # pragma: no cover
    # TODO: switch to argparse.
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--tol", type="float", default=1e-12)
    parser.add_option("--max_iter", type="int", default=1000000)
    parser.add_option("--cutoff", type="float", default=0.05)
    parser.add_option("--pstats", type="string", default=None)
    parser.add_option("--float_format", type="string", default="10.8")
    opts, args = parser.parse_args(parseargs)
    return opts, args


if __name__ == "__main__":  # pragma: no cover
    opts, args = parsing(sys.argv[1:])
    wham = WHAM(tol=opts.tol, max_iter=opts.max_iter, cutoff=opts.cutoff)
    wham.float_format = opts.float_format
    df = wham.load_files(args)

    # keep around some stuff to allow us to do benchmarking and profiling
    if opts.pstats is not None:
        import cProfile
        import time
        start = time.time()
        cProfile.run("print wham.wham_bam_histogram(df)", opts.pstats)
        print(time.time() - start)
    else:
        wham_hist = wham.wham_bam_histogram(df)
        print(wham_hist.to_string(
            header=False,
            float_format=lambda x: "{:10.8f}".format(x)
        ))
