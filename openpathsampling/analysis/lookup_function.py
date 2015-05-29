import pandas as pd
import numpy as np
class LookupFunction(object):
    """
    Interpolation between datapoints.

    Parameters
    ----------
    ordinate : iterable of numbers
        values for the ordinate
    abscissa : iterable of numbers
        values for the abscissa

    Iteration and numpy ufuncs work on the values. Callable with any number.

    Notes
    -----
        Largely, this class mimics an immutable dictionary, except instead
        of implementing __getitem__, we use the __call__ function. If you
        call a number that is in the dictionary, you get exactly that
        number. If you call a number that it not in the dictionary, the get
        the linear interpolation/extrapolation for that number based on the
        dictionary values.
    """
    def __init__(self, ordinate, abscissa):
        self.pairs = { }
        for (x,y) in zip(ordinate, abscissa):
            self.pairs[x] = y
        self.sorted_ordinates = np.array(sorted(self.pairs.keys()))
        self._values = np.array([self.pairs[x] for x in self.sorted_ordinates])

    def keys(self):
        """
        Return the (ordered) list of ordinates
        """
        return list(self.sorted_ordinates)

    def values(self):
        """
        Return the list of values (ordered by ordinate)
        """
        return self._values

    @property
    def x(self):
        """
        Property to return the ordinates
        """
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
        """Return a pandas.Series representation of data points"""
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

        while (i < nvals and xvals[i] < value):
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
    """
    Simple mean and std for a group of LookupFunctions.

    The mean and std from this are, themselves, LookupFunctions, and so can
    interpolate between included values. Calling the group acts as calling
    the mean. __getitem__, __setitem__, and append act on the list of
    functions.

    Parameters
    ----------
    functions : list of LookupFunctions
        the functions included
    use_x : "shared" (default), "all", or list of numbers
        the values to consider as the ordinates. If "shared", includes only
        values which appear in all the functions. If "all", includes all
        values which appear in any function. A list of numbers will use that
        list as the ordinate values.

    Notes
    -----
        The choice of `use_x` is very important for the calculation of the
        mean and standard deviation: if you use "shared", then you only
        calculate the mean/std at points where all functions have measured
        values. If you use "all", you will include points which are
        interpolated/extrapolated, instead of measured. In the current
        implementation, there is no way to get a mean/std with different
        numbers of contributions at each point, depending on whether the
        point has a measurement or is an extrapolation.
    """
    def __init__(self, functions, use_x="shared"):
        self.functions = functions
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
        """Standard deviation."""
        std = []
        for val in self.x:
            std.append(
                np.array([fcn(val) for fcn in self.functions]).std()
            )
        return LookupFunction(self.x, std)

    @property
    def mean(self):
        """Mean."""
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
