"""
Tools for resampling functions that output pandas.DataFrame objects.

Typically, you use these tools in 3 steps:

1. Create resampling groups of data, using, e.g., BlockResampling
2. Create a function that maps a list of input data into the desired output
   DataFrame
3. Create a ResamplingStatistics object using the `input` from step 1 and
   the `function` from step 2.
"""

import numpy as np
import itertools
import pandas as pd

import logging
logger = logging.getLogger(__name__)

# NOTE: there may be a better way to do this, by converting the results to
# numpy arrays and using the numpy functions. However, you'd have to be very
# careful that the rows and columns still correspond to the same things,
# i.e., that there's no permutation of index order. Using pandas protects us
# from such problems. Alternatively, it may be possible to bunch these into
# a pandas.Panel, and use pandas functions. All requires some performance
# testing; the advantange of the approach below is that when we want both
# mean and std, the mean can be calculated first and passed to the std to
# speed it up. (However, performance probably won't matter much for the
# things we'll be doing).

def mean_df(objects):
    """Basic calculation of mean (average) of a list of DataFrames.

    Parameters
    ----------
    objects : list of pandas.DataFrame
        the DataFrames to calculate the mean (NB: technically, these don't
        have to be DataFrames. They must be closed on the `sum` operation
        and over division by a float.

    Returns
    -------
    pandas.DataFrame :
        the mean of each element in the DataFrame
    """
    return sum(objects) / float(len(objects))

def std_df(objects, mean_x=None):
    """Standard deviation of a list of DataFrames.

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
    n_obj = float(len(objects))
    sq = [o**2 for o in objects]
    variance = mean_df(sq) - mean_x**2
    return variance.applymap(np.sqrt)

class ResamplingStatistics(object):
    """
    Contains and organizes resampled statistics.

    Attributes
    ----------
    results : list of pandas.DataFrame
        the result DataFrames after applying the function
    mean : pandas.DataFrame
        DataFrame containing mean of the result DataFrames
    std : pandas.DataFrame
        DataFrame containing standard deviation of the result DataFrames

    Parameters
    ----------
    function : callable
        the function to apply the statistics to; must take one item from
        the list `inputs` and return a pandas.DataFrame
    inputs : list
        each element of inputs is can be used as input to `function`
    """
    def __init__(self, function, inputs):
        self.function = function
        self.inputs = inputs
        self.results = [self.function(inp) for inp in self.inputs]
        self._mean = None
        self._std = None
        self._sorted_series = None

    @property
    def mean(self):
        if self._mean is None:
            self._mean  = mean_df(self.results)
        return self._mean

    @property
    def std(self):
        if self._std is None:
            self._std = std_df(self.results, mean_x=self.mean)
        return self._std

    @property
    def index(self):
        return self.mean.index

    @property
    def columns(self):
        return self.mean.columns

    @property
    def sorted_series(self):
        if self._sorted_series is None:
            self._sorted_series = {
                loc: pd.Series(df.loc[loc]
                               for df in self.results).sort_values()
                for loc in itertools.product(self.index, self.columns)
            }
        return self._sorted_series

    def percentile(self, percent):
        """Percentile, using Nearest Rank method.

        Parameters
        ----------
        percent : float
            the percentile desired

        Returns
        -------
        pd.DataFrame
            the DataFrame nearest to that percentile
        """
        n_entries = len(self.results)
        rank = min(int(percent / 100.0 * n_entries), n_entries - 1)
        df = pd.DataFrame(index=self.index, columns=self.columns)
        for idx in self.index:
            for col in self.columns:
                df.loc[idx, col] = self.sorted_series[(idx, col)].iloc[rank]
        return df

class BlockResampling(object):
    """Select samples according to block resampling.

    If neither n_blocks nor n_per_block are set (as is the default
    behavior) then n_blocks=20 is used. The blocks are always of the same
    size, if the number of samples doesn't divide evenly, then the extra
    samples are placed in the `unassigned` attribute.

    Parameters
    ----------
    all_samples : list
        list of all samples
    n_blocks : int
        number of blocks (resampling sets)
    n_per_block : int
        number of samples per block
    """
    def __init__(self, all_samples, n_blocks=None, n_per_block=None):
        self.n_total_samples = len(all_samples)
        if n_blocks is None and n_per_block is None:
            n_blocks = 20
        if n_blocks is None and n_per_block is not None:
            n_blocks = self.n_total_samples // n_per_block
        elif n_blocks is not None and n_per_block is None:
            n_per_block = self.n_total_samples // n_blocks

        self.n_blocks = n_blocks
        self.n_per_block = n_per_block
        self.blocks = [all_samples[i*n_per_block:(i+1)*n_per_block]
                       for i in range(n_blocks)]
        self.unassigned = all_samples[n_blocks*n_per_block:]
        self.n_resampled = self.n_total_samples - len(self.unassigned)
