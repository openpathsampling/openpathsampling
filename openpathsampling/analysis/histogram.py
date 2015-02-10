import numpy as np

class Histogram(object):
    """Wrapper for numpy.histogram with additional conveniences.

    In addition to the behavior in numpy.histogram, this provides a few
    additional calculations, as well as behavior that allows for better
    interactive use (tricks to assist cachingby libraries using it, etc.)
    """
    def __init__(self, n_bins=40, bin_width=None, bin_range=None):
        """Creates the parameters for the histogram.

        Either `n_bins` or `bin_width` must be given. If both are given,
        `bin_width` overrides `n_bins`.
        """
        pass

    def add_data_to_histogram(self, data, weights):
        """Adds more data to an existing histogram"""
        pass

    def build_from_data(self, data, weights):
        """Builds the histogram based on `data`
        
        Note
        ----
        Repeating this results in recalculating the histogram. This is the
        expected behavior; in using this, you should check if the histogram
        parameters have changed from a previous run (using
        `compare_parameters`) and you should be aware whether your data has
        changed.
        """
        pass

    def compare_parameters(self, other):
        """Return true if `other` has the same bin parameters as `self`.

        Useful for checking whether a histogram needs to be rebuilt.
        """
        # None returns false: use that as a quick test
        if other == None:
            return False
        match = (
            (self.n_bins == other.n_bins) and
            (self.bin_width == other.bin_width) and 
            (self.bin_range[0] == other.bin_range[0]) and 
            (self.bin_range[1] == other.bin_range[1])
        )
        return match

    def _total_weight(self):
        tot = 0.0
        for bincount in self.histogram:
            tot += bincount
        return tot

    # Yes, the following could be cached. No, I don't think it is worth it.
    # Keep in mind that we need a separate cache for each one that we build,
    # and that typically it will take almost no time to build one of these.
    # Adding caching complicates the code for no real benefit.

    def normalized(self, raw_probability=False):
        """Return normalized version of histogram.

        By default (`raw_probability` false), this returns the histogram
        normalized by its integral (according to rectangle-rule
        integration). If `raw_probability` is true, this returns the
        histogram normalized by the sum of the bin counts, with no
        consideration of the bin widths.
        """
        pass

    def cumulative(self, maximum=None):
        """Cumulative from the left: number of values less than bin value.
        """
        pass

    def reverse_cumulative(self, maximum=None):
        """Cumulative from the right: number of values greater than bin value.
        """
        pass
