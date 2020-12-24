import openpathsampling as paths
import openpathsampling.netcdfplus as netcdfplus
import collections
import copy

from functools import partial

def _cv_max_func(trajectory, cv):
    return max(cv(trajectory))

def _netcdfplus_cv_max(cv):
    return paths.netcdfplus.FunctionPseudoAttribute(
        name="max " + cv.name,
        key_class=paths.Trajectory,
        f=_cv_max_func,
        cv=cv
    ).with_diskcache(allow_incomplete=True)


def _simstore_cv_max(cv):
    from openpathsampling.experimental.simstore import StorableFunction
    return StorableFunction(
        func=_cv_max_func,
        cv=cv
    ).named("max " + cv.name)


class InterfaceSet(netcdfplus.StorableNamedObject):
    """List of volumes representing a set of interfaces, plus metadata.

    Implements (immutable) list API, such that the InterfaceSet can act like
    a list of the interface volumes.

    Parameters
    ----------
    volumes : list of :class:`.Volume`
        volumes representing the interfaces
    cv : :class:`.CollectiveVariable`
        order parameter for this interface set
    lambdas : list
        values associated with the CV at each interface
    cv_max : :class:`.netcdfplus.PseudoAttribute`
        function to calculate the maximum value of the CV (if not given,
        will be generated from ``cv`` if that is given)
    """
    # This is a class variable because more than one interface set may use
    # the same CV (with different interface values) -- common in testing.
    _cv_max_dict = {}
    simstore = False

    def __init__(self, volumes, cv=None, lambdas=None, cv_max=None,
                 direction=None):
        super(InterfaceSet, self).__init__()
        self.volumes = volumes
        self.cv = cv

        known_cvs = list(self._cv_max_dict.keys())
        if cv_max is None:
            if self.cv in known_cvs:
                self.cv_max = self._cv_max_dict[self.cv]
                # NOTE: If KeyError occurs, see InterfaceSet._reset()
            elif self.cv is not None:
                factory = {True: _simstore_cv_max,
                           False: _netcdfplus_cv_max}[self.simstore]
                self.cv_max = factory(cv)
            else:
                self.cv_max = None
        else:
            self.cv_max = cv_max

        if self.cv is not None and self.cv not in known_cvs:
            self._cv_max_dict[self.cv] = self.cv_max

        self.lambdas = lambdas
        self.direction = direction
        if direction is None and lambdas is not None:
            # we guess based on the values of lambda
            # if they aren't monotone, we can't tell and return 0
            count = len(lambdas)-1
            increasing = sum([lambdas[i+1] >= lambdas[i]
                              for i in range(len(lambdas)-1)]) == count
            decreasing = sum([lambdas[i+1] <= lambdas[i]
                              for i in range(len(lambdas)-1)]) == count

            if increasing:
                self.direction = 1
            elif decreasing:
                self.direction = -1
        if self.direction is None:
            self.direction = 0

        self._set_lambda_dict()

    @classmethod
    def _reset(cls):
        """Reset the max_cv_dict.

        Used in testing, when CVs of the same name get created multiple
        times as part of tests. This is only a problem if you have (in
        memory) multiple CV objects with the same name. That can occur in
        tests (where the entire test suite is run on one import) but should
        never happen in normal practice.
        """
        cls._cv_max_dict = {}

    def _set_lambda_dict(self):
        vlambdas = self.lambdas
        if vlambdas is None:
            vlambdas = [None]*len(self.volumes)
        self._lambda_dict = {vol: lmbda
                             for (vol, lmbda) in zip(self.volumes, vlambdas)}

    def get_lambda(self, volume):
        """Lambda (value of the CV) associated with a given interface volume

        Parameters
        ----------
        volume : :class:`.Volume`
            the interface volume

        Returns
        -------
        float or int
            the value of the CV associated with the interface
        """
        return self._lambda_dict[volume]

    def _slice_dict(self, slicer):
        dct = self.to_dict()
        dct['volumes'] = self.volumes[slicer]
        try:
            dct['lambdas'] = self.lambdas[slicer]
        except TypeError:
            dct['lambdas'] = self.lambdas
        return dct

    def __len__(self):
        return len(self.volumes)

    def __getitem__(self, key):
        result = self.volumes[key]
        if type(result) is list:
            return self.from_dict(self._slice_dict(key))
        else:
            return result

    def __iter__(self):
        return iter(self.volumes)

    def __contains__(self, item):
        return item in self.volumes

    def __reversed__(self):
        return self.volumes.__reversed__()


class GenericVolumeInterfaceSet(InterfaceSet):
    """Abstract class for InterfaceSets for CVRange-based volumes.

    Subclasses act as factories for interface volumes, as well as holding
    the metadata about them.

    Parameters
    ----------
    cv : :class:`.CollectiveVariable`
        the collective variable for this interface set
    minvals : float or int or list of float or list of int
        the minimum value(s) for the interface set
    maxvals : float or int or list of float or list of int
        the maximum value(s) for the interface set
    intersect_with : :class:`.Volume`
        output volumes will be intersected (`&`) with this.
    volume_func : callable, returns :class:.`Volume`, takes minval, maxval
        the function to create the interface volume based on the CV.
        Typically the differentiating factor of subclasses.
    """
    def __init__(self, cv, minvals, maxvals, intersect_with, volume_func):
        if intersect_with is None:
            intersect_with = paths.FullVolume()
        self.intersect_with = intersect_with
        self.minvals = minvals
        self.maxvals = maxvals

        minvs, maxvs, direction = self._sanitize_input(minvals, maxvals)
        lambdas = {1: maxvs, -1: minvs, 0: None}[direction]
        volumes = [self.intersect_with & volume_func(minv, maxv)
                   for (minv, maxv) in zip(minvs, maxvs)]
        super(GenericVolumeInterfaceSet, self).__init__(
            volumes=volumes, cv=cv, lambdas=lambdas, direction=direction)
        self._set_volume_func(volume_func)

    def _slice_dict(self, slicer):
        dct = super(GenericVolumeInterfaceSet, self)._slice_dict(slicer)
        try:
            dct['minvals'] = self.minvals[slicer]
        except TypeError:
            dct['minvals'] = self.minvals
        try:
            dct['maxvals'] = self.maxvals[slicer]
        except TypeError:
            dct['maxvals'] = self.maxvals
        return dct

    def _set_volume_func(self, volume_func):
        if self.direction == 0:
            self.volume_func = volume_func
        elif self.direction > 0:
            self.volume_func = lambda maxv: volume_func(self.minvals, maxv)
        elif self.direction < 0:
            self.volume_func = lambda minv: volume_func(minv, self.maxvals)

    def to_dict(self):
        dct = {'cv': self.cv,
               'minvals': self.minvals,
               'maxvals': self.maxvals,
               'intersect_with': self.intersect_with,
               'lambdas': self.lambdas,
               'direction': self.direction,
               'volumes': self.volumes}
        try:
            dct.update({'cv_max': self.cv_max})
        except AttributeError:  # pragma: no cover
            # reverse compatibility; deprecated in 0.9.5, remove in 2.0
            pass

        return dct

    def _load_from_dict(self, dct):
        self.cv = dct['cv']
        try:
            self.cv_max = dct['cv_max']
        except KeyError:  # pragma: no cover
            # reverse compatibility; deprecated in 0.9.4, remove in 2.0
            pass
        self.minvals = dct['minvals']
        self.maxvals = dct['maxvals']
        self.intersect_with = dct['intersect_with']
        self.lambdas = dct['lambdas']
        self.direction = dct['direction']
        self.volumes = dct['volumes']
        self._set_lambda_dict()

    @staticmethod
    def _sanitize_input(minvals, maxvals):
        """Normalizes the input of minvals and maxvals.

        Parameters
        ----------
        minvals : float or int or list of float or list of int
            the minimum value(s) for the interface set
        maxvals : float or int or list of float or list of int
            the maximum value(s) for the interface set

        Returns
        -------
        minvals : list
            minimum values as a list
        maxvals : list
            maximum values as a list
        direction : 1, -1, or 0
            whether the maximum value are increasing (1), the minimum values
            are decreasing (-1), or it is unclear (0). "Unclear" can happen
            if both are changing or if neither are changing.
        """
        direction = 0
        try:
            len_min = len(minvals)
        except TypeError:
            len_min = 1
            minvals = [minvals]
        try:
            len_max = len(maxvals)
        except TypeError:
            len_max = 1
            maxvals = [maxvals]
        if len_min == len_max:
            # check if all elements of each list matches its first element
            if minvals.count(minvals[0]) == len_min:
                direction += 1
            if maxvals.count(maxvals[0]) == len_max:
                direction += -1
            # this approach means that if multiple vals are equal (for some
            # drunken reason, you decided to have a bunch of equivalent
            # volumes?) we return that we can't tell the direction
        elif len_max > len_min == 1:
            direction = 1
        elif len_min > len_max == 1:
            direction = -1
        else:
            raise RuntimeError("Can't reconcile array lengths: "
                               + str(minvals) + ", " + str(maxvals))

        minvs = minvals
        maxvs = maxvals
        if len_min == 1:
            minvs = minvs*len(maxvs)
        if len_max == 1:
            maxvs = maxvs*len(minvs)
        return minvs, maxvs, direction

    def new_interface(self, lambda_i):
        """Creates a new interface at lambda_i.

        Note
        ----
        This only returns the interface; it does *not* add it to this
        interface set.

        Parameters
        ----------
        lambda_i : float or int
            the value of the CV to associated with the new interface.

        Returns
        -------
        :class:`.Volume`
            new interface volume

        Raises
        ------
        TypeError
            If the volume_func requires both a minimum and a maximum value
            of lambda. Message will say "<lambda>() takes exactly 2
            arguments (1 given)".
        """
        return self.intersect_with & self.volume_func(lambda_i)



class VolumeInterfaceSet(GenericVolumeInterfaceSet):
    """InterfaceSet based on CVDefinedVolume.

    Parameters
    ----------
    cv : :class:`.CollectiveVariable`
        the collective variable for this interface set
    minvals : float or int or list of float or list of int
        the minimum value(s) for the interface set
    maxvals : float or int or list of float or list of int
        the maximum value(s) for the interface set
    intersect_with : :class:`.Volume`
        output volumes will be intersected (`&`) with this.
    """
    def __init__(self, cv, minvals, maxvals, intersect_with=None):
        volume_func = lambda minv, maxv: paths.CVDefinedVolume(cv, minv, maxv)
        super(VolumeInterfaceSet, self).__init__(cv, minvals, maxvals,
                                                 intersect_with,
                                                 volume_func)

    @staticmethod
    def from_dict(dct):
        interface_set = VolumeInterfaceSet.__new__(VolumeInterfaceSet)
        interface_set._load_from_dict(dct)
        volume_func = lambda minv, maxv: paths.CVDefinedVolume(
            interface_set.cv, minv, maxv
        )
        super(InterfaceSet, interface_set).__init__()
        interface_set._set_volume_func(volume_func)
        return interface_set


class PeriodicVolumeInterfaceSet(GenericVolumeInterfaceSet):
    """InterfaceSet based on PeriodicCVDefinedVolume.

    Parameters
    ----------
    cv : :class:`.CollectiveVariable`
        the collective variable for this interface set
    minvals : float or int or list of float or list of int
        the minimum value(s) for the interface set
    maxvals : float or int or list of float or list of int
        the maximum value(s) for the interface set
    period_min : float (optional)
        minimum of the periodic domain
    period_max : float (optional)
        maximum of the periodic domain
    intersect_with : :class:`.Volume`
        output volumes will be intersected (`&`) with this.
    """
    def __init__(self, cv, minvals, maxvals, period_min=None,
                 period_max=None, intersect_with=None):
        volume_func = lambda minv, maxv: paths.PeriodicCVDefinedVolume(
            cv, minv, maxv, period_min, period_max
        )
        self.period_min = period_min
        self.period_max = period_max
        super(PeriodicVolumeInterfaceSet, self).__init__(cv, minvals,
                                                         maxvals,
                                                         intersect_with,
                                                         volume_func)

    def to_dict(self):
        dct = super(PeriodicVolumeInterfaceSet, self).to_dict()
        dct['period_min'] = self.period_min
        dct['period_max'] = self.period_max
        return dct

    @classmethod
    def from_dict(cls, dct):
        interface_set = cls.__new__(cls)

        interface_set._load_from_dict(dct)
        interface_set.period_min = dct['period_min']
        interface_set.period_max = dct['period_max']
        volume_func = lambda minv, maxv: paths.PeriodicCVDefinedVolume(
            interface_set.cv,
            minv, maxv,
            interface_set.period_min,
            interface_set.period_max
        )
        super(InterfaceSet, interface_set).__init__()
        interface_set._set_volume_func(volume_func)
        return interface_set
