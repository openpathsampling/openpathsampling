import openpathsampling as paths
import openpathsampling.netcdfplus as netcdfplus

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
    """
    def __init__(self, volumes, cv=None, lambdas=None):
        self.volumes = volumes
        self.cv = cv
        self.lambdas = lambdas
        try:
            self.direction = lambdas[-1] >= lambdas[0]
        except TypeError:
            self.direction = 0

        vlambdas = lambdas
        if vlambdas is None:
            vlambdas = [None]*len(volumes)
        self._lambda_dict = {vol: lmbda 
                             for (vol, lmbda) in zip(volumes, vlambdas)}

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

    def __len__(self):
        return len(self.volumes)

    def __getitem__(self, key):
        return self.volumes[key]

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

        minvs, maxvs, direction = self._sanitize_input(minvals, maxvals)
        lambdas = {1: maxvs, -1: minvs, 0: None}[direction]
        volumes = [self.intersect_with & volume_func(minv, maxv)
                   for (minv, maxv) in (minvs, maxvs)]
        super(self, GenericVolumeInterfaceSet).__init__(volumes, cv, lambdas)

        if direction == 0:
            self.volume_func = volume_func
        elif direction > 0:
            self.volume_func = lambda maxv : volume_func(minvals, maxv)
        elif direction < 0:
            self.volume_func = lambda minv : volume_func(minv, maxvals)

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
        """
        return self.intersect_with & self.volume_func(lambda_i)


class VolumeInterfaceSet(GenericVolumeInterfaceSet):
    """InterfaceSet based on CVRangeVolume.

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
        volume_func = lambda minv, maxv: paths.CVRangeVolume(cv, minv, maxv)
        super(self, VolumeInterfaceSet).__init__(cv, minvals, maxvals,
                                                 intersect_with,
                                                 volume_func)


class PeriodicVolumeInterfaceSet(GenericVolumeInterfaceSet):
    """InterfaceSet based on CVRangeVolumePeriodic.

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
        volume_func = lambda minv, maxv: paths.CVRangeVolumePeriodic(
            cv, minv, maxv, period_min, period_max
        )
        super(self, VolumeInterfaceSet).__init__(cv, minvals, maxvals,
                                                 intersect_with,
                                                 volume_func)


