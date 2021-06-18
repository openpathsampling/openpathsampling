'''
Created on 03.09.2014

@author: Jan-Hendrik Prinz, David W.H. Swenson
'''

from . import range_logic
import abc
from openpathsampling.netcdfplus import StorableNamedObject
import numpy as np
import warnings

# TODO: Make Full and Empty be Singletons to avoid storing them several times!

def join_volumes(volume_list, name=None):
    """
    Make the union of a list of volumes. (Useful shortcut.)

    Parameters
    ----------
    volume_list : list of :class:`openpathsampling.Volume`
        the list to be joined together
    name : str or callable
        string for name, or callable that creates string for name from
        ``volume_list``

    Returns
    -------
    :class:`openpathsampling.UnionVolume`
        the union of the elements of the list, or EmptyVolume if list is
        empty
    """
    volume = EmptyVolume()
    # EmptyVolume is smart and knows its OR just takes the other
    for vol in volume_list:
        volume = volume | vol
    if name is not None:
        try:
            name_str = name(volume_list)
        except TypeError:
            name_str = name
        volume = volume.named(name_str)
    return volume


class Volume(StorableNamedObject):
    """
    A Volume describes a set of snapshots
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self):
        super(Volume, self).__init__()

    @abc.abstractmethod
    def __call__(self, snapshot):
        '''
        Returns `True` if the given snapshot is part of the defined Region
        '''
        return False # pragma: no cover

    def __str__(self):
        '''
        Returns a string representation of the volume
        '''
        return 'volume' # pragma: no cover

    __hash__ = StorableNamedObject.__hash__

    def __or__(self, other):
        if self is other:
            return self
        elif type(other) is EmptyVolume:
            return self
        elif type(other) is FullVolume:
            return other
        else:
            return UnionVolume(self, other)

    def __xor__(self, other):
        if self is other:
            return EmptyVolume()
        elif type(other) is EmptyVolume:
            return self
        elif type(other) is FullVolume:
            return ~ self
        else:
            return SymmetricDifferenceVolume(self, other)

    def __and__(self, other):
        if self is other:
            return self
        elif type(other) is EmptyVolume:
            return other
        elif type(other) is FullVolume:
            return self
        else:
            return IntersectionVolume(self, other)

    def __sub__(self, other):
        if self is other:
            return EmptyVolume()
        elif type(other) is EmptyVolume:
            return self
        elif type(other) is FullVolume:
            return EmptyVolume()
        else:
            return RelativeComplementVolume(self, other)

    def __invert__(self):
        return NegatedVolume(self)

    def __eq__(self, other):
        return str(self) == str(other)

    def __ne__(self, other):
        return not self == other


class VolumeCombination(Volume):
    """
    Logical combination of volumes.

    This should be treated as an abstract class. For storage purposes, use
    specific subclasses in practice.
    """
    def __init__(self, volume1, volume2, fnc, str_fnc):
        super(VolumeCombination, self).__init__()
        self.volume1 = volume1
        self.volume2 = volume2
        self.fnc = fnc
        self.sfnc = str_fnc

    def __call__(self, snapshot):
        # short circuit following JHP's implementation in ensemble.py
        a = self.volume1(snapshot)
        res_true = self.fnc(a, True)
        res_false = self.fnc(a, False)
        if res_false == res_true:
            return res_true
        else:
            b = self.volume2(snapshot)
            return self.fnc(a, b)
        #return self.fnc(self.volume1.__call__(snapshot),
                        #self.volume2.__call__(snapshot))

    def __str__(self):
        return '(' + self.sfnc.format(str(self.volume1), str(self.volume2)) + ')'

    def to_dict(self):
        return {'volume1': self.volume1, 'volume2': self.volume2}


class UnionVolume(VolumeCombination):
    """ "Or" combination (union) of two volumes."""
    def __init__(self, volume1, volume2):
        super(UnionVolume, self).__init__(
            volume1=volume1,
            volume2=volume2,
            fnc=lambda a, b: a or b,
            str_fnc='{0} or {1}'
        )


class IntersectionVolume(VolumeCombination):
    """ "And" combination (intersection) of two volumes."""
    def __init__(self, volume1, volume2):
        super(IntersectionVolume, self).__init__(
            volume1=volume1,
            volume2=volume2,
            fnc=lambda a, b: a and b,
            str_fnc='{0} and {1}'
        )


class SymmetricDifferenceVolume(VolumeCombination):
    """ "Xor" combination of two volumes."""
    def __init__(self, volume1, volume2):
        super(SymmetricDifferenceVolume, self).__init__(
            volume1=volume1,
            volume2=volume2,
            fnc=lambda a, b: a ^ b,
            str_fnc='{0} xor {1}'
        )


class RelativeComplementVolume(VolumeCombination):
    """ "Subtraction" combination (relative complement) of two volumes."""
    def __init__(self, volume1, volume2):
        super(RelativeComplementVolume, self).__init__(
            volume1=volume1,
            volume2=volume2,
            fnc=lambda a, b: a and not b,
            str_fnc='{0} and not {1}'
        )


class NegatedVolume(Volume):
    """Negation (logical not) of a volume."""
    def __init__(self, volume):
        super(NegatedVolume, self).__init__()
        self.volume = volume

    def __call__(self, snapshot):
        return not self.volume(snapshot)

    def __str__(self):
        return '(not ' + str(self.volume) + ')'


class EmptyVolume(Volume):
    """Empty volume: no snapshot can satisfy"""
    def __init__(self):
        super(EmptyVolume, self).__init__()

    def __call__(self, snapshot):
        return False

    def __and__(self, other):
        return self

    def __or__(self, other):
        return other

    def __xor__(self, other):
        return other

    def __sub__(self, other):
        return self

    def __invert__(self):
        return FullVolume()

    def __str__(self):
        return 'empty'


class FullVolume(Volume):
    """Volume which all snapshots can satisfy."""
    def __init__(self):
        super(FullVolume, self).__init__()

    def __call__(self, snapshot):
        return True

    def __invert__(self):
        return EmptyVolume()

    def __and__(self, other):
        return other

    def __or__(self, other):
        return self

    def __xor__(self, other):
        return ~ other

    def __sub__(self, other):
        return ~ other

    def __str__(self):
        return 'all'


class CVDefinedVolume(Volume):
    """
    Volume defined by a range of a collective variable `collectivevariable`.

    Contains all snapshots `snap` for which `lamba_min <=
    collectivevariable(snap)` and `lambda_max > collectivevariable(snap)`.

    Parameters
    ----------
    collectivevariable : :class:`.CollectiveVariable`
        the CV to base the volume on
    lambda_min : float
        minimum value of the CV
    lambda_max : float
        maximum value of the CV
    """
    def __init__(self, collectivevariable, lambda_min=0.0, lambda_max=1.0):
        super(CVDefinedVolume, self).__init__()
        self.collectivevariable = collectivevariable
        try:
            self.lambda_min = lambda_min.__float__()
        except AttributeError:
            self.lambda_min = float(lambda_min)

        try:
            self.lambda_max = lambda_max.__float__()
        except AttributeError:
            self.lambda_max = float(lambda_max)

        self._cv_returns_iterable = None  # used to raise warnings

    # Typically, the logical combinations are only done once. Because of
    # this, it is worth passing these through a check to speed up the logic.

    # To get all the usefulness of the range logic in a subclass, all you
    # should need to override is _copy_with_new_range (so that it inits any
    # extra info the subclass carries) and range_and/or/sub, so that they
    # return the correct behavior for the new subclass. Everything else
    # comes for free.

    @property
    def default_name(self):
        return (str(self.lambda_min) + "<"
                + str(self.collectivevariable.name) + "<"
                + str(self.lambda_max))

    def _copy_with_new_range(self, lmin, lmax):
        """Shortcut to make a CVDefinedVolume with all parameters the same as
        this one except the range. This is useful for the range logic when
        dealing with subclasses: just override this function to copy extra
        information.
        """
        return CVDefinedVolume(self.collectivevariable, lmin, lmax)

    @staticmethod
    def range_and(amin, amax, bmin, bmax):
        return range_logic.range_and(amin, amax, bmin, bmax)
    @staticmethod
    def range_or(amin, amax, bmin, bmax):
        return range_logic.range_or(amin, amax, bmin, bmax)
    @staticmethod
    def range_sub(amin, amax, bmin, bmax):
        return range_logic.range_sub(amin, amax, bmin, bmax)

    def _lrange_to_Volume(self, lrange):
        """Takes results from one of the range_logic functions and returns
        the appropriate Volume.

        Parameters
        ----------
        lrange : None or 1 or list of 2-tuples
            Key to the volume to be returned: None returns the EmptyVolume, 1
            returns self, and a list of 2-tuples is __or__'d as (min,max) to
            make a VolumeCombinations

        Returns
        -------
        Volume
            appriate volume according to lrange

        Raises
        ------
        ValueError
            if the input lrange is not an allowed value
        """
        if lrange is None:
            return EmptyVolume()
        elif lrange == 1:
            return self
        elif lrange == -1:
            return FullVolume()
        elif len(lrange) == 1:
            return self._copy_with_new_range(lrange[0][0], lrange[0][1])
        elif len(lrange) == 2:
            return UnionVolume(
                self._copy_with_new_range(lrange[0][0], lrange[0][1]),
                self._copy_with_new_range(lrange[1][0], lrange[1][1])
            )
        else:
            raise ValueError(
                "lrange value not understood: {0}".format(lrange)
            )  # pragma: no cover

    def __and__(self, other):
        if (type(other) is type(self) and
                self.collectivevariable == other.collectivevariable):
            lminmax = self.range_and(self.lambda_min, self.lambda_max,
                                other.lambda_min, other.lambda_max)
            return self._lrange_to_Volume(lminmax)
        else:
            return super(CVDefinedVolume, self).__and__(other)

    def __or__(self, other):
        if (type(other) is type(self) and
                self.collectivevariable == other.collectivevariable):
            lminmax = self.range_or(self.lambda_min, self.lambda_max,
                               other.lambda_min, other.lambda_max)
            return self._lrange_to_Volume(lminmax)
        else:
            return super(CVDefinedVolume, self).__or__(other)

    def __xor__(self, other):
        if (type(other) is type(self) and
                self.collectivevariable == other.collectivevariable):
            # taking the shortcut here
            return (self | other) - (self & other)
        else:
            return super(CVDefinedVolume, self).__xor__(other)

    def __sub__(self, other):
        if (type(other) is type(self) and
                self.collectivevariable == other.collectivevariable):
            lminmax = self.range_sub(self.lambda_min, self.lambda_max,
                            other.lambda_min, other.lambda_max)
            return self._lrange_to_Volume(lminmax)
        else:
            return super(CVDefinedVolume, self).__sub__(other)

    def _is_iterable(self, val):
        try:
            # simtk.Quantity erroneously allows iter, so use len
            # besides, CVs shouldn't return generators
            _ = len(val)
        except TypeError:
            return False
        else:
            cv = self.collectivevariable
            warnings.warn("The CV '" + str(cv.name) + "' returns an "
                          "iterable. This may lead to problem in analysis.")
            return True

    def _get_cv_float(self, snapshot):
        val = self.collectivevariable(snapshot)
        if self._cv_returns_iterable is None:
            self._cv_returns_iterable = self._is_iterable(val)
        return val.__float__()

    def __call__(self, snapshot):
        l = self._get_cv_float(snapshot)

        # we explicitly test for infinity to allow the user to
        # define `lambda_min/max='inf'` also when using units
        # a simtk unit cannot be compared to a python infinite float
        if self.lambda_min != float('-inf') and self.lambda_min > l:
            return False

        if self.lambda_min != float('inf') and self.lambda_max <= l:
            return False

        return True

    def __str__(self):
        return '{{x|{2}(x) in [{0:g}, {1:g}]}}'.format(
            self.lambda_min, self.lambda_max, self.collectivevariable.name)

class PeriodicCVDefinedVolume(CVDefinedVolume):
    """
    As with `CVDefinedVolume`, but for a periodic order parameter.

    Defines a Volume containing all states where collectivevariable, a periodic
    function wrapping into the range [period_min, period_max], is in the
    given range [lambda_min, lambda_max].

    Attributes
    ----------
    period_min : float (optional)
        minimum of the periodic domain
    period_max : float (optional)
        maximum of the periodic domain
    """

    _excluded_attr = ['wrap']

    def __init__(
            self, collectivevariable, lambda_min=0.0, lambda_max=1.0,
            period_min=None, period_max=None):
        super(PeriodicCVDefinedVolume, self).__init__(collectivevariable,
                                                    lambda_min, lambda_max)
        self.period_min = period_min
        self.period_max = period_max
        if (period_min is not None) and (period_max is not None):
            self._period_shift = period_min
            self._period_len = period_max - period_min
            if self.lambda_max - self.lambda_min > self._period_len:
                raise Exception("Range of volume larger than periodic bounds.")
            elif self.lambda_max-self.lambda_min == self._period_len:
                # this is only the case that we really have a FullVolume
                self.lambda_min = period_min
                self.lambda_max = period_max
                # hack: better to create factory, returning FullVolume
                # this hack: https://stackoverflow.com/questions/38541015/
                class MonkeyPatch(type(self)):
                    def __call__(self, *arg, **kwarg):
                        return True
                self.__class__ = MonkeyPatch
            else:
                self.lambda_min = self.do_wrap(lambda_min)
                self.lambda_max = self.do_wrap(lambda_max)
            self.wrap = True
        else:
            self.wrap = False

    def do_wrap(self, value):
        """Wraps `value` into the periodic domain."""

        # this looks strange and mimics the modulo operation `%` while
        # being fully compatible for simtk numbers and plain python as well
        # working for ints and floats.
        val = value - self._period_shift

        # little trick to check for positivity without knowing the the units
        # or if it actually has units

        if val > val * 0:
            return value - int(val / self._period_len) * self._period_len
        else:
            wrapped = value + int((self._period_len - val) / self._period_len) \
                * self._period_len
            if wrapped >= self._period_len:
                wrapped -= self._period_len

            return wrapped

    # next few functions add support for range logic
    def _copy_with_new_range(self, lmin, lmax):
        return PeriodicCVDefinedVolume(self.collectivevariable, lmin, lmax,
                                    self.period_min, self.period_max)

    @staticmethod
    def range_and(amin, amax, bmin, bmax):
        return range_logic.periodic_range_and(amin, amax, bmin, bmax)

    @staticmethod
    def range_or(amin, amax, bmin, bmax):
        return range_logic.periodic_range_or(amin, amax, bmin, bmax)

    @staticmethod
    def range_sub(amin, amax, bmin, bmax):
        return range_logic.periodic_range_sub(amin, amax, bmin, bmax)

    def __invert__(self):
        # consists of swapping max and min
        return PeriodicCVDefinedVolume(self.collectivevariable,
                                    self.lambda_max, self.lambda_min,
                                    self.period_min, self.period_max
                                   )

    def __call__(self, snapshot):
        l = self._get_cv_float(snapshot)
        if self.wrap:
            l = self.do_wrap(l)
        if self.lambda_min > self.lambda_max:
            return l >= self.lambda_min or l < self.lambda_max
        else:
            return self.lambda_min <= l < self.lambda_max

    def __str__(self):
        if self.wrap:
            fcn = 'x|({0}(x) - {2:g}) % {1:g} + {2:g}'.format(
                        self.collectivevariable.name,
                        self._period_len, self._period_shift)
            if self.lambda_min < self.lambda_max:
                domain = '[{0:g}, {1:g}]'.format(
                        self.lambda_min, self.lambda_max)
            else:
                domain = '[{0:g}, {1:g}] union [{2:g}, {3:g}]'.format(
                        self._period_shift, self.lambda_max,
                        self.lambda_min, self._period_shift+self._period_len)
            return '{'+fcn+' in '+domain+'}'
        else:
            return '{{x|{2}(x) [periodic] in [{0:g}, {1:g}]}}'.format(
                        self.lambda_min, self.lambda_max,
                        self.collectivevariable.name)


class VoronoiVolume(Volume):
    '''
    Volume given by a Voronoi cell specified by a set of centers

    Parameters
    ----------
    collectivevariable : MultiRMSDCV
        must be an MultiRMSDCV collectivevariable that returns several RMSDs
    state : int
        the index of the center for the chosen voronoi cell

    Attributes
    ----------
    collectivevariable : collectivevariable
        the collectivevariable object
    state : int
        the index of the center for the chosen voronoi cell

    '''

    def __init__(self, collectivevariable, state):
        super(VoronoiVolume, self).__init__()
        self.collectivevariable = collectivevariable
        self.state = state

    def cell(self, snapshot):
        '''
        Returns the index of the voronoicell snapshot is in

        Parameters
        ----------
        snapshot : :class:`opensampling.engines.BaseSnapshot`
            the snapshot to be tested

        Returns
        -------
        int
            index of the voronoi cell
        '''
        distances = self.collectivevariable(snapshot)
        min_val = 1000000000.0
        min_idx = -1
        for idx, d in enumerate(distances):
            if d < min_val:
                min_val = d
                min_idx = idx

        return min_idx

    def __call__(self, snapshot, state=None):
        '''
        Returns `True` if snapshot belongs to voronoi cell in state

        Parameters
        ----------
        snapshot : :class:`opensampling.engines.BaseSnapshot`
            snapshot to be tested
        state : int or None
            index of the cell to be tested. If `None` (Default) then the
            internal self.state is used

        Returns
        -------
        bool
            returns `True` is snapshot is on the specified voronoi cell

        '''
        if state is None:
            state = self.state

        return self.cell(snapshot) == state


# class VolumeFactory(object):
    # @staticmethod
    # def _check_minmax(minvals, maxvals):
        # # if one is an integer, convert it to a list
        # if type(minvals) == int or type(minvals) == float:
            # if type(maxvals) == list:
                # minvals = [minvals]*len(maxvals)
            # else:
                # raise ValueError("minvals is a scalar; maxvals is not a list")
        # elif type(maxvals) == int or type(maxvals) == float:
            # if type(minvals) == list:
                # maxvals = [maxvals]*len(minvals)
            # else:
                # raise ValueError("maxvals is a scalar; minvals is not a list")

        # if len(minvals) != len(maxvals):
            # raise ValueError("len(minvals) != len(maxvals)")
        # return (minvals, maxvals)

    # @staticmethod
    # def CVRangeVolumeSet(op, minvals, maxvals):
        # # TODO: clean up to only use min_i or max_i in name if necessary
        # minvals, maxvals = VolumeFactory._check_minmax(minvals, maxvals)
        # myset = []
        # for (min_i, max_i) in zip(minvals, maxvals):
            # volume = CVDefinedVolume(op, min_i, max_i)
            # myset.append(volume)
        # return myset

    # @staticmethod
    # def CVRangeVolumePeriodicSet(op, minvals, maxvals,
                                # period_min=None, period_max=None):
        # minvals, maxvals = VolumeFactory._check_minmax(minvals, maxvals)
        # myset = []
        # for i in range(len(maxvals)):
            # myset.append(PeriodicCVDefinedVolume(op, minvals[i], maxvals[i],
                                              # period_min, period_max))
        # return myset
