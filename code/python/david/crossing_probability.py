#!/usr/bin/env python
import os, sys
from histogram import read_file
from series_stats import add_series_to_set

def parsing(sysargs):
    import optparse
    parser = optparse.OptionParser()
    return parser.parse_args(sysargs)


def series_set_crossing_lambdas(series_set):
    lambda_i = { }
    prev_set = series_set[min(series_set.keys())]
    prev_lambda = 0
    for lmbda in sorted(series_set.keys()):
        for i in range(len(series_set[lmbda])):
            if (prev_set[i]==1.0) and (series_set[lmbda][i]!=1.0):
                lambda_i[i] = prev_lambda
        prev_lambda = lmbda
        prev_set = series_set[lmbda]
    lambda_i[len(lambda_i)] = max(series_set.keys())
    return lambda_i
    


def series_set_crossing_prob(series_set):
    crossing_probs = { }
    lambda_i = series_set_lambdas(series_set)

    print "#", lambda_i
    print "# 1.0",
    for i in sorted(lambda_i.keys()):
        if i<max(lambda_i.keys()):
            print series_set[lambda_i[i+1]][i],
    print

    scaling = 1.0
    i=0

    # this is the version that segments the sampling: between each pair of
    # interfaces, only the data from the previous interface set is used
    for lmbda in sorted(series_set.keys()):
        if (lmbda==lambda_i[i+1] and lmbda!=max(series_set.keys())):
            scaling *= series_set[lambda_i[i+1]][i]
            i += 1 
        crossing_probs[lmbda] = scaling*series_set[lmbda][i]
        print lmbda, crossing_probs[lmbda]
    
    # other version
    #crossing_probs2 = { }
    #for lmbda in sorted(series_set.keys()):
    #    crossing_probs2[lmbda] = 0.0
    #    for i in range(len(series_set[lmbda])):
    #        scaling = 1.0
    #        if (lmbda < lambda_i[i]):
    #            scaling = 0.0
    #       else:
    #           j=i
    #           while (j>0):
    #               scaling *= series_set[lambda_i[j+1]][i]
    #               j -= 1
    #       crossing_probs2[lmbda] += scaling * series_set[lmbda][i]
        #print lmbda, crossing_probs2[lmbda]

    #for lmbda in sorted(series_set.keys()):
    #    print lmbda, crossing_probs[lmbda], crossing_probs2[lmbda]
        
        
    return crossing_probs

if __name__ == "__main__":
    opt, args = parsing(sys.argv)
    all_series = { }
    for fname in args[1:]:
        keys, vals = read_file(fname)
        add_series_to_set(keys, vals, all_series)

    crossing_probs = series_set_crossing_prob(all_series)
    
    #for key in crossing_probs.keys().sort():
    #    print key, crossing_probs[key]
