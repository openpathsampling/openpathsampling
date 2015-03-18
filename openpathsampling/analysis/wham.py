#!/usr/bin/env python

# wham.py
# The Weighted Histogram Analysis Method (WHAM) for combining histograms
# from several files.
import math
import pandas as pd

import logging
logger = logging.getLogger(__name__)
#from read_acf import read_acf
#from series_stats import add_series_to_set
#from crossing_probability import series_set_crossing_lambdas

class WHAM(object):
    """
    Weighted Histogram Analysis Method
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
        """
        self.tol = tol
        self.max_iter = max_iter
        self.cutoff = cutoff
        self.hists = { }
        self.lnZ = []
        self.sum_hist = []
        self.nt = []
        self.nhists = 0
        self.keys = []

    # loads files into the dictionary
    def load_files(self,fnames):
        """Load a file or files into the internal structures.

        Requires either pandas or something else with pandas-like read_table
        and concat functions.
        """
        frames = []
        try:
            for fname in fnames:
                frames.append(pd.read_table(fname, index_col=0, sep=" ",
                                            usecols=[0,1]))
        except TypeError:
            frame.append(pd.read_table(fnames, index_col=0, sep=" ",
                                       usecols=[0,1]))
        df = pd.concat(frames, axis=1)
        self.load_from_dataframe(df)
        #xvals, yvals = read_acf(fname)
        #add_series_to_set(xvals, yvals, self.hists)
        #self.nhists += 1

    def load_from_dataframe(self, df):
        """Loads from a pandas-compatible data frame.

        The data frame does not need to be from pandas, but needs to quack
        like it is. In particular, it needs a .index property which gives
        the list of row keys, and a .loc[key] property which gives the row
        associated with a given key.
        """
        self.keys = list(df.index)
        for key in self.keys:
            self.hists[key] = list(df.loc[key])
            if len(self.hists[key]) != self.nhists:
                if self.nhists != 0:
                    raise Warning("Inequal number of histograms for WHAM")
                else:
                    self.nhists = len(self.hists[key])



    # modifies the dictionary to ignore outside the window
    def prep(self):
        # find the max; calculate nt
        self.keys = sorted(self.hists.keys())
        for i in range(self.nhists):
            max_i=0
            self.nt.append(0.0)
            for key in self.keys:
                if (self.hists[key][i] > max_i):
                    max_i = self.hists[key][i]

            for k in range(len(self.keys)):
                key = self.keys[k]
                if (self.hists[key][i] < self.cutoff * max_i):
                    self.hists[key][i] = 0
                self.nt[i] += self.hists[key][i]

        # calculate sum_hist
        kill_keys = []
        for k in range(len(self.keys)):
            key = self.keys[k]
            self.sum_hist.append(0.0)
            self.sum_hist[k] = math.fsum(self.hists[key])
            if (self.sum_hist[k] == 0.0):
                kill_keys.append(key)
        
        # remove keys with sum_hist==0
        for deadkey in kill_keys:
            self.sum_hist.remove(0.0)
            self.keys.remove(deadkey)

        #print self.keys
        #print self.sum_hist
        #print self.nt
        #print math.fsum(self.nt), math.fsum(self.sum_hist)

        return

    def clean_leading_ones(self):
        """
        Removes leading ones from input cumulative histograms.
        """
        keys = sorted(self.hists.keys())
        for s in range(self.nhists):
            i=1 # ignore the line of zeros before 
            i=0 # new version doesn't have extra line of zeros
            while (     (self.hists[keys[i]][s]==1.0) \
                    and (self.hists[keys[i+1]][s]==1.0) \
                    and (i < len(keys)-1) ):
                self.hists[keys[i]][s] = 0.0
                i += 1
        return

    def guess_lnZ(self):
        """
        Prepares a guess of ln(Z) based on crossing probabilities.

        A trivial guess should also work, but this isn't very expensive.
        """
        # taken from
        # $DYNQ/scripts/crossing_probability.py:series_set_crossing_lambdas
        lambda_i = { }
        prev_set = self.hists[min(self.hists.keys())]
        prev_lambda = 0
        for lmbda in sorted(self.hists.keys()):
            for i in range(len(self.hists[lmbda])):
                if (prev_set[i]==1.0) and (self.hists[lmbda][i]!=1.0):
                    lambda_i[i] = prev_lambda
            prev_lambda = lmbda
            prev_set = self.hists[lmbda]
        lambda_i[len(lambda_i)] = max(self.hists.keys())

        #lambda_i = series_set_crossing_lambdas(self.hists)
        scaling = 1.0
        crossing_prob = []
        for key in sorted(lambda_i.keys()):
            #print key, lambda_i[key]
            if (key > 1):
                crossing_prob.append(self.hists[lambda_i[key-1]][key-2])
            elif (key > 0):
                crossing_prob.append(1.0)

        scaled = []
        scaling = 1.0
        for val in crossing_prob:
            self.lnZ.append(math.log(scaling*val))
            scaled.append(scaling*val)
            scaling *= val

        return

    # wham iterations; returns the WHAM lnZ weights
    def generate_lnZ(self):
        diff = self.tol+1
        iteration = 0
        while (diff > self.tol) and (iteration < self.max_iter):
            lnZnew = []
            for histnum in range(self.nhists):
                newZi=0
                for k in range(len(self.keys)):
                    key = self.keys[k]
                    if (self.hists[key][histnum] > 0):
                        mysum = 0
                        for j in range(self.nhists):
                            weight = 0.0
                            if (self.hists[key][j] > 0):
                                weight = 1.0
                            mysum += weight*self.nt[j]/math.exp(self.lnZ[j])
                        newZi += self.sum_hist[k] / mysum

                lnZnew.append(math.log(newZi))

            # get error
            diff=0
            for histnum in range(self.nhists):
                diff += abs(self.lnZ[histnum] - lnZnew[histnum])
                self.lnZ[histnum] = lnZnew[histnum] - lnZnew[0]

            iteration += 1

            # check status (mainly for debugging)
            sampling = 1 + self.max_iter # DEBUG
            if (iteration % sampling == 0):
                logger.debug("niteration = " + str(iteration))
                logger.debug("  diff = " + str(diff))
                logger.debug("   lnZ = " + str(self.lnZ))
                logger.debug("lnZnew = " + str(lnZnew))

        logger.info("iterations=" + str(iteration) + "diff=" + str(diff))
        logger.info("       lnZ=" + str(self.lnZ))
        self.convergence = (iteration, diff)


    def wham_histogram(self):

        final_hist ={}
        for key in self.keys:
            final_hist[key] = []
            for histnum in range(self.nhists):
                orig = self.hists[key][histnum]
                final_hist[key].append( orig*math.exp(self.lnZ[histnum]) / \
                        self.nt[histnum])

        # generate (unnormalized) final hist
        wham_hist = {}
        for k in range(len(self.keys)):
            key = self.keys[k]
            mysum=0
            for histnum in range(self.nhists):
                weight = 0
                if (self.hists[key][histnum] > 0):
                    weight = 1
                mysum += weight*self.nt[histnum] / math.exp(self.lnZ[histnum])
            if (mysum > 0):
                wham_hist[key] = self.sum_hist[k] / mysum
            else:
                wham_hist[key] = 0 # raise a red flag

        # normalize hist
        norm=0
        for key in self.keys:
            if (wham_hist[key] > norm):
                norm = wham_hist[key]

        for key in self.keys:
            wham_hist[key] /= norm

        return wham_hist

    def wham_bam_histogram(self):
        self.guess_lnZ()
        self.prep()
        self.generate_lnZ()
        hist = self.wham_histogram()
        return hist



def parsing(parseargs):
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--tol", type="float", default=1e-12)
    parser.add_option("--max_iter", type="int", default=1000000)
    parser.add_option("--cutoff", type="float", default=0.05)
    opts, args = parser.parse_args(parseargs)
    return opts, args

def print_dict(adict):
    keys = adict.keys()
    keys.sort()
    for key in keys:
        print key, adict[key]

import sys, os
if __name__ == "__main__":
    opts, args = parsing(sys.argv[1:])
    wham = WHAM(tol=opts.tol, max_iter=opts.max_iter, cutoff=opts.cutoff)
    wham.load_files(args)

    wham.clean_leading_ones()
    wham_hist = wham.wham_bam_histogram()

    print_dict(wham_hist)


