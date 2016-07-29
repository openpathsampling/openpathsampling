import openpathsampling as paths

class InterfaceSet(paths.StorableNamedObject):
    def __init__(self, volumes, cv=None, lambdas=None, direction=1):
        self.volumes = volumes
        self.cv = cv
        self.lambdas = lambdas
        self.direction = direction

        vlambdas = lambdas
        if vlambdas is None:
            vlambdas = [None]*len(volumes)
        self._lambda_dict = {vol: lmbda 
                             for (vol, lmbda) in zip(volumes, vlambdas)}

    def get_lambda(self, vol):
        return self._lambda_dict[vol]

    # TODO: add various magics (contains, iter, etc) so that this acts like
    # a list over its volume objects. Append will have to be different
    # (requires lambda as well as a volume).


class GenericVolumeInterfaceSet(InterfaceSet):
    def __init__(self, cv, minvals, maxvals, intersect_with, volume_func):
        if intersect_with is None:
            intersect_with = paths.FullVolume()
        self.intersect_with = intersect_with

        minvs, maxvs = self._prep_minvals_maxvals(minvals, maxvals)
        direction = self._determine_direction(minvals, maxvals)
        lambdas = {1: maxvs, -1: minvs, 0: None}[direction]
        volumes = [self.intersect_with & volume_func(minv, maxv)
                   for (minv, maxv) in (minvs, maxvs)]

        if direction == 0:
            self.volume_func = volume_func
        elif direction > 0:
            self.volume_func = lambda maxv : volume_func(minvals, maxv)
        elif direction < 0:
            self.volume_func = lambda minv : volume_func(minv, maxvals)

    @staticmethod
    def _determine_direction(minvals, maxvals):
        try:
            len_min = len(minvals)
        except TypeError:
            len_min = 1
        try:
            len_max = len(maxvals)
        except TypeError:
            len_max = 1
        if len_min == len_max:
            return 0
        elif len_max > len_min:
            return 1
        else:
            return -1

    @staticmethod
    def _prep_minvals_maxvals(minvals, maxvals, direction):
        minvs = minvals
        maxvs = maxvals
        if direction > 0:
            minvs = [minvs]*len(maxvs)
        elif direction < 0:
            maxvs = [maxvs]*len(maxvs)
        return minvs, maxvs

    def new_interface(self, lambda_i):
        return self.intersect_with & self.volume_func(lambda_i)


class VolumeInterfaceSet(GenericVolumeInterfaceSet):
    def __init__(self, cv, minvals, maxvals, intersect_with=None):
        volume_func = lambda minv, maxv: paths.CVRangeVolume(cv, minv, maxv)
        super(self, VolumeInterfaceSet).__init__(cv, minvals, maxvals,
                                                 intersect_with,
                                                 volume_func)


class PeriodicVolumeInterfaceSet(GenericVolumeInterfaceSet):
    def __init__(self, cv, minvals, maxvals, period_min=None,
                 period_max=None, intersect_with=None):
        volume_func = lambda minv, maxv: paths.CVRangeVolumePeriodic(
            cv, minv, maxv, period_min, period_max
        )
        super(self, VolumeInterfaceSet).__init__(cv, minvals, maxvals,
                                                 intersect_with,
                                                 volume_func)


