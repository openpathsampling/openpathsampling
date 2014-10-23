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
