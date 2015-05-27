import pandas as pd
import numpy as np
class LookupFunction(object):
    def __init__(self, ordinate, abscissa):
        self.pairs = { }
        for (x,y) in zip(ordinate, abscissa):
            self.pairs[x] = y
        self.sorted_ordinates = np.array(sorted(self.pairs.keys()))
        self._values = np.array([self.pairs[x] for x in self.sorted_ordinates])

    def keys(self):
        return list(self.sorted_ordinates)

    def values(self):
        return self._values

    @property
    def x(self):
        return self.sorted_ordinates


    def __len__(self):
        return len(self.sorted_ordinates)

    def __iter__(self):
        for val in self.values():
            yield val

    # TODO: may need better array behaviors
    def __array__(self, result=None):
        return np.array(self.values())

    def __array_wrap__(self, result, context=None):
        res_arr = np.ndarray.__array_wrap__(self._values, result, context)
        return LookupFunction(self.sorted_ordinates, res_arr)

    def __array_prepare__(self, result, context=None):
        return result

    def series(self):
        # TODO: temp hack until I can get matplotlib to plot natively
        ser = pd.Series(self.values(), self.keys())
        return ser

    def __call__(self, value):
        # only a 1D implementation so far
        i=0
        xvals = self.sorted_ordinates
        nvals = len(xvals)
        if value < xvals[i]:
            # extrapolation TODO: add log warning
            x1 = xvals[0]
            x2 = xvals[1]

        while (xvals[i] < value and i < nvals):
            i += 1

        if i == nvals:
            # extrapolation TODO: add log warning
            x1 = xvals[-2]
            x2 = xvals[-1]
        else:
            # interpolation
            x1 = xvals[i-1]
            x2 = xvals[i]

        y1 = self.pairs[x1]
        y2 = self.pairs[x2]

        y = float(value - x1) / (x2 - x1) * (y2-y1) + y1
        return y



class LookupFunctionGroup(LookupFunction):
    def __init__(self, lookup_functions, use_x="shared"):
        self.functions = lookup_functions
        self.shared_x = set(self.functions[0].x)
        self.all_x = set(self.functions[0].x)
        for fcn in self.functions:
            self.shared_x = self.shared_x & set(fcn.x)
            self.all_x = self.all_x | set(fcn.x)
        
        self.use_x = use_x


    @property
    def use_x(self):
        return self._use_x

    @use_x.setter
    def use_x(self, use_x):
        self._use_x = use_x
        if use_x == "all":
            self.sorted_ordinates = self.all_x
        elif use_x == "shared":
            self.sorted_ordinates = self.shared_x
        else:
            self.sorted_ordinates = use_x

    @property
    def std(self):
        std = []
        for val in self.x:
            std.append(
                np.array([fcn(val) for fcn in self.functions]).std()
            )
        return LookupFunction(self.x, std)

    @property
    def mean(self):
        mean = []
        for val in self.x:
            mean.append(
                np.array([fcn(val) for fcn in self.functions]).mean()
            )
        return LookupFunction(self.x, mean)

    def __call__(self, value):
        return  self.mean(value)


    def __getitem__(self, item):
        return self.functions[item]

    def __setitem__(self, item, value):
        self.functions[item] = value

    def __contains__(self, item):
        return item in self.functions
    
    def append(self, item):
        self.functions.append(item)
