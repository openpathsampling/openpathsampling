'''
@author: David W.H. Swenson

This file consists of the logic functions to combine two sets defined as
ranges of numbers. The functions are such that "and" is the intersection of
two sets, "or" is the union of two sets, and "sub" means that A - B is the
relative complement of B in A (usually denoted A \ B).

The sets are defined by their minima and maxima; for set A this is amin and
amax, and for set B this is bmin and bmax. Since everything in this code is
done in terms of amin, amax, bmin, and bmax, we often just refer to these
points as "a", "A", "b", "B" respectively, using the lower-case for the
minima and upper-case for maxima.

Defining the logic here means that we can explicitly use this logic to
simplify combinations at either the volume level (when the two volumes have
the same order parameter) or at the ensemble level (e.g., for
`LengthEnsemble`s).

All of these functions return results which can then be translated into a
volume according to the rules of the CVDefinedVolume._lrange_to_Volume
function:

          Return Value  | Volume generated
------------------------+--------------------------------------------------
                     1  | `self` for volume A
                    -1  | `FullVolume`
                  None  | `EmptyVolume`
 list with one 2-tuple  | `CVDefinedVolume` with range given by 2-tuple
list with two 2-tuples  | `VolumeCombination` of or'd together
                        | `CVDefinedVolume`s with ranges given by the 2-tuples
------------------------+--------------------------------------------------


PERIODIC RANGES
===============

In the nonperiodic case, the fact that a<A and b<B simplifies the test logic
substantially. However, in the periodic case, this is not necessarily true;
A<a simply means that the volume wraps around the periodic edge, which is
quite possible. This makes the periodic case much more complex.  Another
part of the challenge is that we don't require the user to define the
periodic domain. Because of this, we don't have a trivial way to unwrap the
periodic volumes.

The approach we will use attempts to reduce the periodic problem to the
nonperiodic case whenever possible. The basic outline of it is to order
the `a`, `A`, `b`, and `B` variables. Then, instead of using the actual
values associated with those, we (when possible) call the nonperiodic range
logic using dummy variables from the ordering. So the order a<A<b<B calls
the range logic with amin=0, amax=1, bmin=2, bmax=3, wheras the order
b<a<B<A calls the range logic with amin=1, amax=3, bmin=0, bmax=2. Then we
recover the actual numbers from the range that comes out of that.

In a little more detail:

The first step is to find the sequence of the minima and maxima. In the
periodic case, this sequence is only unique up to the set of cyclic
permutations. For these four elements, we define a canonical cyclic
permutation such that we always have a<A, and if possible, we have b<B. If
both of these are satified, we prefer a<b. The algorithm to generate this
works by first assuming that `a` is less than all the others. This gives us
3*2*1 permutations, which we categorize into the following groups:

(1) Satisfy all rules: aAbB abAB abBA
(2) breaks b<B rule, but `b` is largest: aABb, aBAb
(3) breaks b<B rule, but `b` is not largest: aBbA

Once a<A and b<B, we can simply use the logic from the nonperiodic case.
Group (1) satisfies those rules trivially. Group (2) can be made to satisfy
those rules through a cyclic permutation, i.e., aABb=>baAB and aBAb=>baBA,
and then that can be send to the nonperiodic functions. Only group (3),
which only contains one element, needs to be treated as a special case for
the periodic range logic functions.

The function `periodic_ordering` returns the order such that `order[0]` is
the order index for `amin`, `order[1]` for `amax`, etc., with the order
index being integers between 0 and 3 inclusive.

Then, if possible, we call the nonperiodic range logic. This returns the
appropriate combination. For example, if we called `and` on a function which
gave the order abAB=(amin=0, amax=2, bmin=1, bmax=3), `range_and` woudl
return [(1,2)], saying that the correct range is between 1 and 2, i.e,
between bmin and amax.

We put this, along the with `order` we got from periodic ordering and a
dictionary of order-to-value, into the function `recover_periodic_range`. It
first uses `order` on the tuple to switch the order numbers used in the
range logic back to labels corresponding to amin, amax, etc. Then the
dictionary maps those labels to the actual input values, and out comes the
correct range tuple.

Total no-brainer, right?
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
    amin : float
        minimum of first range
    amax : float
        maximum of first range
    bmin : float
        minimum of second range
    bmax : float
        maximum of second range

    Returns
    -------
    list of int 0-3
        Order index of amin, amax, bmin, bmax in that order; i.e. the return
        value [0, 2, 1, 3] means amin < bmin < amax < bmax; amin in order
        spot 0, amax in 2, bmin in 1, bmax in 3.
    """
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
        out = [out[3]] + out[slice(0, 3)]
    # at this point we have a canonically ordered list of the letter a, A,
    # b, and B.
    final = [out.index(a) for a in ['a', 'A', 'b', 'B']]
    return final

def recover_periodic_range(lrange, order, adict, aa_is_full=False):
    # aa_is_full: [a:a] represents FullEnsemble
    if lrange == None or lrange == 1:
        return lrange
    else:
        retval = []
        for pair in lrange:
            opair = [order.index(pval) for pval in pair]
            opair = [order.index(pval) for pval in pair]
            mytup = tuple(map(adict.get, opair))
            if (mytup[0] == mytup[1]):
                # periodic ordering sometimes leads to [a:a] ranges; for
                # and/sub, this should be the empty volume; for or it should
                # be the full ensemble
                if aa_is_full == True:
                    return -1
                # otherwise, this is the empty ensemble, so ignore
            else:
                retval.append(mytup)
        if len(retval) == 0:
            # I don't think we can get here. If we do, correct behavior is:
            return None # pragma: no cover
        else:
            return retval

def periodic_range_and(amin, amax, bmin, bmax):
    adict = {0: amin, 1: amax, 2: bmin, 3: bmax}
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
        return recover_periodic_range(or_res, order, adict, aa_is_full=True)

def periodic_range_sub(amin, amax, bmin, bmax):
    adict = {0 : amin, 1 : amax, 2 : bmin, 3 : bmax}
    if amin == bmin and amax == bmax:
        return None
    elif amin == bmax and bmin == amax:
        return 1
    order = periodic_ordering(amin, amax, bmin, bmax)
    if order == [0, 3, 2, 1]: # aBbA
        special_res = [(1, 2)] # order[1] to order[2]
        return recover_periodic_range(special_res, order, adict)
    else:
        sub_res = range_sub(order[0], order[1], order[2], order[3])
        return recover_periodic_range(sub_res, order, adict)

