import pandas as pd
import numpy as np
class LookupFunction(object):
    def __init__(self, ordinate, abscissa):
        self.pairs = { }
        for (x,y) in zip(ordinate, abscissa):
            self.pairs[x] = y
        self.sorted_ordinates = np.array(sorted(self.pairs.keys()))

    def keys(self):
        return list(self.sorted_ordinates)

    def values(self):
        return np.array([self.pairs[x] for x in self.sorted_ordinates])

    @property
    def x(self):
        return self.sorted_ordinates


    def __len__(self):
        return len(self.sorted_ordinates)

    # TODO: may need better array behaviors
    def __array__(self, result=None):
        return np.array(self.values())

    def __array_wrap__(self, result, context=None):
        return LookupFunction(self.values(), result)

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
    def __init__(self, lookup_functions):
        self.functions = lookup_functions
        pass

    def std(self):
        pass

    def __getitem__(self):
        pass

    def __setitem__(self):
        pass
        
