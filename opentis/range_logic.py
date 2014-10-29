'''
@author: David W.H. Swenson
'''

def range_and(amin, amax, bmin, bmax):
    if amin == bmin and amax == bmax:
        return 1
    if amin < bmin:
        lmin = bmin
    else:
        lmin = amin
    if amax > bmax:
        lmax = bmax
    else:
        lmax = amax
    if lmin > lmax:
        return None
    else:
        return [(lmin, lmax)]

def range_or(amin, amax, bmin, bmax):
    if amin == bmin and amax == bmax:
        return 1
    if amin < bmin:
        lmin = amin
        gmin = bmin
    else:
        lmin = bmin
        gmin = amin
    if amax > bmax:
        lmax = amax
        gmax = bmax
    else:
        lmax = bmax
        gmax = amax
    if gmax > gmin:
        return [(lmin, lmax)]
    else:
        return [(amin, amax), (bmin, bmax)]


def range_sub(amin, amax, bmin, bmax):
    """A - B"""
    if amin == bmin and amax == bmax:
        return None
    if bmax < amin or bmin > amax:
        return 1
    if amin < bmin:
        if amax > bmax:
            return [(amin, bmin), (bmax, amax)]
        else:
            return [(amin, bmin)]
    else:
        if amax < bmax:
            return None
        else:
            return [(bmax, amax)]

# I make no claim that this is the fastest algorithm to solve this problem.
# However, it is simple, and it is designed to allow me to re-use the code
# above for periodic systems as well. This is unlikely to be a performance
# bottleneck, so readability and testability rate much higher than speed.
def periodic_ordering(amin, amax, bmin, bmax):
    """Figures out the order of the permutation that maps the minima and
    maxima to their order, in canonical form (amin<amax, bmin<bmax if
    possible).

    Parameters
    ----------

    Return
    -------
        : list of int 0-3
        Order index of amin, amax, bmin, bmax in that order; i.e. the return
        value [0, 2, 1, 3] means amin < bmin < amax < bmax; amin in order
        spot 0, amax in 2, bmin in 1, bmax in 3. 
    """
    dict1 = {amin : 'a', amax : 'A', bmin : 'b', bmax : 'B'}
    dict2 = {'a' : amin, 'A' : amax, 'b' : bmin, 'B' : bmax}
    order = ['a']
    # put the labels in the increasing order, starting at amin
    for label in ('A', 'b', 'B'):
        i = 0
        val = dict2[label]
        while i < len(order):
            if val < dict2[order[i]]:
                order.insert(i, label)
                break
            i += 1
        if label not in order:
            order.append(label)
    # Canonical order is 'a' always before 'A', and if possible, 'b' before
    # 'B', and 'a' before 'B' (in that order of priority). This defines a
    # unique member within the set of cyclic permutations; find that cyclic
    # permutation.
    idx0 = order.index('a')
    out = []
    for i in range(4):
        out.append(order[(idx0+i) % 4])
    if out[3] == 'b':
        out = [out[3]] + out[slice(0,3)]
    final = [out.index(a) for a in ['a','A','b','B']]
    return final

def recover_periodic_range(lrange, order, adict):
    if lrange == None or lrange == 1:
        return lrange
    else:
        retval = []
        for pair in lrange:
            opair = [ order.index(pval) for pval in pair ]
            retval.append(tuple(map(adict.get, opair)))
        return retval

def periodic_range_and(amin, amax, bmin, bmax):
    adict = {0 : amin, 1 : amax, 2 : bmin, 3 : bmax}
    if amin == bmin and amax == bmax:
        return 1
    if amin == bmax and bmin == amax:
        return None
    order = periodic_ordering(amin, amax, bmin, bmax)
    if order == [0, 3, 2, 1]: # aBbA
        special_res = [(0, 1), (2, 3)]
        return recover_periodic_range(special_res, order, adict)
    else:
        and_res = range_and(order[0], order[1], order[2], order[3])
        return recover_periodic_range(and_res, order, adict)

def periodic_range_or(amin, amax, bmin, bmax):
    adict = {0 : amin, 1 : amax, 2 : bmin, 3 : bmax}
    if amin == bmin and amax == bmax:
        return 1
    elif amin == bmax and bmin == amax:
        return -1
    order = periodic_ordering(amin, amax, bmin, bmax)
    if order == [0, 3, 2, 1]: # aBbA
        return -1
    else:
        or_res = range_or(order[0], order[1], order[2], order[3])
        return recover_periodic_range(or_res, order, adict)

def periodic_range_sub(amin, amax, bmin, bmax):
    adict = {0 : amin, 1 : amax, 2 : bmin, 3 : bmax}
    if amin == bmin and amax == bmax:
        return None
    elif amin == bmax and bmin == amax:
        return None
    order = periodic_ordering(amin, amax, bmin, bmax)
    if order == [0, 3, 2, 1]: # aBbA
        special_res = [(3, 2)]
        return recover_periodic_range(special_res, order, adict)
    else:
        sub_res = range_sub(order[0], order[1], order[2], order[3])
        return recover_periodic_range(sub_res, order, adict)

