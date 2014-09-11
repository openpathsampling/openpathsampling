#!/usr/bin/env python
import os, sys

class histogram(object):
    
    def __init__(self):
        self.nbins = 0
        self.min = 0
        self.max = 0
        self.tol = 1e-10
        self.keys = []
        self.vals = []
        self.tot_val = 0.0      # for the mean
        self.tot_valsq = 0.0    # for the stddev
        self.count = 0.0
        return
    #TODO


tol = 1e-10

def add_to_hist(hist, val, weight, max, min, nbins):
    bin_num = 1
    bin_max = min
    bin_delta = (max-min)/nbins
    while ((val >= bin_max) and (bin_num < nbins+2)): # +2 for out of range
        bin_num = bin_num + 1
        bin_max = bin_max + bin_delta

    hist[bin_num-1] = hist[bin_num-1] + weight
    #if (bin_num == nbins+1):
        #print val, weight
    return

def find_bin(hist, val, weight, max, min, nbins):
    bin_num = 1
    bin_max = min  # first box is everything before 
    bin_delta = (max-min)/nbins
    # the weird way of saying val > bin_max catches a really weird error
    # that seems to be related to numbers that are very much almost equal
    # but have differences on the order of 10^-17 -- should also put exactly
    # equal into the bin that starts with that value
    while ( (val >= bin_max*(1-tol) ) and (bin_num < nbins+2) ):
        bin_num += 1
        bin_max += bin_delta

    return bin_num

# NOTE: this works to add to bins with SMALLER values than the value being
# tested -- normally I would think of cumulative histograms as adding to
# bins with LARGER values (so from a Gaussian, this should give us erf(-x)
# instead of erf(x))
def add_to_cumhist(hist, val, weight, max, min, nbins):
    bin_num = 1
    bin_max = min
    bin_delta = (max-min)/nbins
    while ((val >= bin_max*(1-tol)) and (bin_num < nbins+2)):
        bin_num = bin_num + 1
        bin_max = bin_max + bin_delta
        hist[bin_num-1] = hist[bin_num-1] + weight

    return

def read_next_line(ffile):
    from re import sub, split, search
    line = ffile.readline()
    line = sub('\#.*', '', line)
    line = sub('^\s*', '', line)
    line = sub('\s*$', '', line)
    splitter = split('\s+', line)
    return splitter


def read_file(filename):
    from re import sub, split, search
    myfile = sys.stdin if filename=="" else open(filename, "r")
    vals = []
    weights = []
    
    line = read_next_line(myfile)
    lnum = 0
    while ((len(line) > 0) and (line[0] != '')):  
        vals.append(float(line[0]))
        weights.append(float(line[1]))
        line = read_next_line(myfile)

    #for line in myfile:
    #    line = sub('\#.*', '', line)
    #    line = sub('^\s*', '', line)
    #    line = sub('\s*$', '', line)
    #    splitter = split('\s+', line)
    #    if (len(splitter) > 1):
    #        vals.append(float(splitter[0]))
    #        weights.append(float(splitter[1]))

    myfile.close()

    return vals, weights

def parsing(sysargs):
    import optparse
    parser=optparse.OptionParser()
    parser.add_option('--min', type='float')
    parser.add_option('--max', type='float')
    parser.add_option('--nbins', type='int')
    parser.add_option('-P', '--probabilities', action="store_true",
                      dest="probabilities", help="modifies normalization")
    parser.add_option('-C', '--cumulative', action="store_true",
                      dest="cumulative")
    parser.add_option('-N', '--normalize', action="store_true", 
                      dest="normalize", help=
"""Normalize: if cumulative option is on, normalized to the maximum value.
If probabilities option is on, normalize such that the sum of the histogram
points is 1. Otherwise, normalize such that the integral of the histogram is
1.  """)
    return parser.parse_args(sysargs)


if __name__ == "__main__":
    opt, args = parsing(sys.argv)
    fname = "" if len(args)==1 else args[1]

    # build an empty histogram
    hist = []
    for i in range(opt.nbins+2):
        hist.append(0.0)

    # test with, e.g.,
    #print hist


    myfile = sys.stdin if fname=="" else open(fname, "r")
    line = read_next_line(myfile)
    while ((len(line) > 0) and (line[0]!='')):
        val = float(line[0])
        weight = float(line[1])
        if (opt.cumulative):
            add_to_cumhist(hist,val,weight, opt.max,opt.min,opt.nbins)
        else:
            add_to_hist(hist,val,weight, opt.max,opt.min,opt.nbins)
        line = read_next_line(myfile)

            


    #(vals, weights) = read_file(fname)
    #if (opt.cumulative):
    #    for i in range(len(vals)):
    #        add_to_cumhist(hist,vals[i],weights[i],opt.max,opt.min,opt.nbins)
    #else:
    #    for i in range(len(vals)):
    #        add_to_hist(hist,vals[i],weights[i],opt.max,opt.min,opt.nbins)

    sum_weights = 0.0
    for i in range(len(hist)):
        sum_weights = sum_weights + hist[i]

    norm = 1.0
    if (opt.normalize):
        if (opt.cumulative):
            # normalize to the maximum value
            norm = 0
            for testval in hist:
                if testval > norm:
                    norm = testval
        elif (opt.probabilities):
            norm = sum_weights 
        else:
            norm = sum_weights * ((opt.max-opt.min) / opt.nbins)
        if norm==0.0:
            print "Empty histogram"
            sys.exit()

    for i in range(len(hist)):
        print opt.min+(i-1)*(opt.max-opt.min)/opt.nbins, hist[i] / norm

    #print "#sum weights =", sum_weights
