"""
Tools for resampling functions that output pandas.DataFrame objects.
"""

def mean_df(objects):
    """Basic calculation of mean (average) of a list of objects.

    Parameters
    ----------
    objects : list of objects
        the objects to be averaged; must be summable (i.e., sum(objects)
        must work and return an object of the same type) and dividable by
        float (i.e., dividing by a float must work and return an object of
        the same type).
    """
    return sum(objects) / float(len(objects))

def std_df(objects, mean_x=None):
    """
    Parameters
    ----------
    objects : list of pandas.DataFrame
        the DataFrames to calculate the standard deviation
    mean_x : pandas.DataFrame
        (optional) pre-calculated mean of the `objects` list. If None
        (default), then the mean will be calculated.

    Returns
    -------
    pandas.DataFrame
        the standard deviation of each element in the dataframe
    """
    if mean_x is None:
        mean_x = mean_df(objects)
    sq = [o**2 for o in objects]
    variance = mean_df(sq) - mean_x**2
    return variance.applymap(math.sqrt)

class ResamplingStatistics(object):
    """
    Parameters
    ----------
    function : callable
        the function to apply the statistics to.
        Must take one item from the list `inputs` and return a pandas.DataFrame
    inputs : list
        each element of inputs is can be used as input to `function`
    """
    def __init__(self, function, inputs):
        self.function = function
        self.inputs = inputs
        self.results = [self.function(inp) for inp in self.inputs]
        self.mean = mean_df(self.results)
        self.std = std_df(self.results, mean=self.mean)

    def percentile_range(self, min_percentile, max_percentile):
        pass

class BlockResampling(object):
    def __init__(self, all_samples, n_blocks=20, n_samples_per_block=None):
        self.n_samples = len(all_samples)
        if n_samples_per_block is not None:
            n_blocks = self.n_samples / n_samples_per_block
        else:
            n_samples_per_block = self.n_samples / n_blocks

        self.blocks = [all_samples[i:i+n_samples_per_block]
                       for i in range(n_blocks)]
