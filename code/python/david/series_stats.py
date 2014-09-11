#!/usr/bin/python

from math import sqrt

def avg_stddev(series):
    tot = 0.0
    tot_sq = 0.0
    for val in series:
        tot += val
        tot_sq += val*val
    avg = tot / len(series)
    stddev = sqrt(tot_sq/len(series) - avg*avg)
    return (avg, stddev)

def add_series_to_set(time, series, series_set):
    for i in range(len(time)):
        t = time[i]
        if (series_set.has_key(t)):
            series_set[t].append(series[i])
        else:
            series_set[t] = [ series[i] ]

    return series_set

def series_set_stats(series_set):
    #from numpy import array
    avg_series = []
    for key in series_set.keys():
        #serii = array(series_set[key]) # yeah, I know it's not a word
        #series_mean = serii.mean()
        #series_stddev = serii.std()
        series_mean, series_stddev = avg_stddev(series_set[key])
        avg_series.append( (key, series_mean, series_stddev) )

    avg_series.sort()
    return avg_series

def series_set_rmsd(series_set, exact): 
    # note that exact must be a series_set with the series of interest in
    # the first position, and for which the set of keys is a superset of the
    # keys in the series_set argument specified above
    rmsd = { }
    for key in series_set.keys():
        armsd = 0.0
        for val in series_set[key]:
            armsd += (val-exact[key][0])**2
        rmsd[key] = sqrt(armsd / len(series_set[key]))

    return rmsd


def parsing(sysargs):
    import optparse
    parser=optparse.OptionParser()
    return parser.parse_args(sysargs)

import sys
if __name__ == "__main__":
    opts, args = parsing(sys.argv[1:])
    all_series = { }
    from read_acf import read_acf
    for file in args:
        time, series = read_acf(file)
        add_series_to_set(time, series, all_series)

    avg = series_set_stats(all_series)
    for i in range(len(avg)):
        print avg[i][0], avg[i][1], avg[i][2]


