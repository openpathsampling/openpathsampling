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
        return (lmin, lmax)

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
        return (lmin, lmax)
    else:
        return 2


def range_sub(amin, amax, bmin, bmax):
    """A - B"""
    if amin == bmin and amax == bmax:
        return None
    if bmax < amin or bmin > amax:
        return 1
    if amin < bmin:
        if amax > bmax:
            return 2
        else:
            return (amin, bmin)
    else:
        if amax < bmax:
            return None
        else:
            return (bmax, amax)

# I make no claim that this is the fastest algorithm to solve this problem.
# However, it is simple, and it is designed to allow me to re-use the code
# above for periodic systems as well. This is unlikely to be a performance
# bottleneck, so readability and testability rate much higher than speed.
def periodic_ordering(amin, amax, bmin, bmax):
    print
    dict1 = {amin : 'a', amax : 'A', bmin : 'b', bmax : 'B'}
    dict2 = {'a' : amin, 'A' : amax, 'b' : bmin, 'B' : bmax}
    order = ['a']
    # put the labels in the increasing order, starting at amin
    for val in (amax, bmin, bmax):
        i = 0
        label = dict1[val]
        while i < len(order):
            if val < dict2[order[i]]:
                order.insert(i, label)
                break
            i += 1
        if label not in order:
            order.append(label)
    # Canonical order is 'a' before 'A' and 'b' before 'B'; with only these
    # four elements, there is always a cyclic permutation where this is
    # true. Find it:
    idx0 = order.index('a')
    out = []
    for i in range(4):
        out.append(order[(idx0+i) % 4])
    if out[3] == 'b':
        out = [out[3]] + out[slice(0,3)]
    print out
    final = [out.index(a) for a in ['a','A','b','B']]
    print final


