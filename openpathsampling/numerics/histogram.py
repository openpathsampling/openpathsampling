import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from .lookup_function import LookupFunction, VoxelLookupFunction
import collections
import warnings
from functools import reduce


class SparseHistogram(object):
    """
    Base class for sparse-based histograms.

    Parameters
    ----------
    bin_widths : array-like
        bin (voxel) size
    left_bin_edges : array-like
        lesser side of the bin (for each direction)
    """
    def __init__(self, bin_widths, left_bin_edges):
        self.bin_widths = np.array(bin_widths)
        if left_bin_edges is None:
            self.left_bin_edges = None
        else:
            self.left_bin_edges = np.array(left_bin_edges)
        self.count = 0
        self.name = None
        self._histogram = None

    def empty_copy(self):
        """Returns a new histogram with the same bin shape, but empty"""
        return type(self)(self.bin_widths, self.left_bin_edges)

    def histogram(self, data=None, weights=None):
        """Build the histogram.

        Parameters
        ----------
        data : list of list of floats
            input data
        weights : list of floats
            weight for each input data point

        Returns
        -------
        collection.Counter :
            copy of the current counter
        """
        if data is None and self._histogram is None:
            raise RuntimeError("histogram() called without data!")
        elif data is not None:
            self._histogram = collections.Counter({})
            return self.add_data_to_histogram(data, weights)
        else:
            return self._histogram.copy()

    @staticmethod
    def sum_histograms(hists):
        # (w, r) = (hists[0].bin_width, hists[0].bin_range)
        # newhist = Histogram(bin_width=w, bin_range=r)
        newhist = hists[0].empty_copy()
        newhist._histogram = collections.Counter({})

        for hist in hists:
            if not newhist.compare_parameters(hist):
                raise RuntimeError
            newhist.count += hist.count
            newhist._histogram += hist._histogram

        return newhist

    def map_to_float_bins(self, trajectory):
        return (np.asarray(trajectory) - self.left_bin_edges) / self.bin_widths

    def map_to_bins(self, data):
        """
        Parameters
        ----------
        data : np.array
            input data

        Returns
        -------
        tuple:
            the bin that the data represents
        """
        # Reshape data to prevent accidental wrong output
        data = np.asarray(data).reshape(self.left_bin_edges.shape)
        return tuple(np.floor((data - self.left_bin_edges) / self.bin_widths))

    def add_data_to_histogram(self, data, weights=None):
        """Adds data to the internal histogram counter.

        Parameters
        ----------
        data : list or list of list
            input data
        weights : list or None
            weight associated with each datapoint. Default `None` is same
            weights for all

        Returns
        -------
        collections.Counter :
            copy of the current histogram counter
        """
        if self._histogram is None:
            return self.histogram(data, weights)
        if weights is None:
            weights = [1.0]*len(data)

        part_hist = sum((collections.Counter({self.map_to_bins(d): w})
                         for (d, w) in zip(data, weights)),
                        collections.Counter({}))

        self._histogram += part_hist
        self.count += len(data) if weights is None else sum(weights)
        return self._histogram.copy()

    @staticmethod
    def _left_edge_to_bin_edge_type(left_bins, widths, bin_edge_type):
        if bin_edge_type == "l":
            return left_bins
        elif bin_edge_type == "m":
            return left_bins + 0.5 * widths
        elif bin_edge_type == "r":
            return left_bins + widths
        elif bin_edge_type == "p":
            pass  # TODO: patches; give the range
        else:
            raise RuntimeError("Unknown bin edge type: " + str(bin_edge_type))

    def xvals(self, bin_edge_type):
        """Position values for the bin

        Parameters
        ----------
        bin_edge_type : 'l' 'm', 'r', 'p'
            type of values to return; 'l' gives left bin edges, 'r' gives
            right bin edges, 'm' gives midpoint of the bin, and 'p' is not
            implemented, but will give vertices of the patch for the bin

        Returns
        -------
        np.array :
            The values of the bin edges
        """
        int_bins = np.array(self._histogram.keys())
        left_bins = int_bins * self.bin_widths + self.left_bin_edges
        return self._left_edge_to_bin_edge_type(left_bins, self.bin_widths,
                                                bin_edge_type)

    def __call__(self, bin_edge_type="m"):
        return VoxelLookupFunction(left_bin_edges=self.left_bin_edges,
                                   bin_widths=self.bin_widths,
                                   counter=self._histogram)

    def normalized(self, raw_probability=False, bin_edge="m"):
        """
        Callable normalized version of the sparse histogram.

        Parameters
        ----------
        raw_probability : bool
            if True, the voxel size is ignored and the sum of the counts
            adds to one. If False (default), the sum of the counts times the
            voxel volume adds to one.
        bin_edge : string
            not used; here for compatibility with 1D versions

        Returns
        -------
        :class:`.VoxelLookupFunction`
            callable version of the normalized histogram
        """
        voxel_vol = reduce(lambda x, y: x.__mul__(y), self.bin_widths)
        scale = voxel_vol if not raw_probability else 1.0
        norm = 1.0 / (self.count * scale)
        counter = collections.Counter({k: self._histogram[k] * norm
                                       for k in self._histogram.keys()})
        return VoxelLookupFunction(left_bin_edges=self.left_bin_edges,
                                   bin_widths=self.bin_widths,
                                   counter=counter)

    def compare_parameters(self, other):
        """Test whether the other histogram has the same parameters.

        Used to check whether we can simply combine these histograms.

        Parameters
        ----------
        other : :class:`.SparseHistogram`
            histogram to compare with

        Returns
        -------
        bool :
            True if these were set up with equivalent parameters, False
            otherwise
        """
        # None returns false: use that as a quick test
        if other is None:
            return False
        if self.left_bin_edges is None or other.left_bin_edges is None:
            # this is to avoid a numpy warning on the next
            return self.left_bin_edges is other.left_bin_edges
        if self.left_bin_edges != other.left_bin_edges:
            return False
        if self.bin_widths != other.bin_widths:
            return False
        return True


class Histogram(SparseHistogram):
    """Wrapper for numpy.histogram with additional conveniences.

    In addition to the behavior in numpy.histogram, this provides a few
    additional calculations, as well as behavior that allows for better
    interactive use (tricks to assist caching by libraries using it, etc.)
    """
    def __init__(self, n_bins=None, bin_width=None, bin_range=None):
        """Creates the parameters for the histogram.

        Either `n_bins` or `bin_width` must be given. If `bin_width` is
        used, then `bin_range` is required. If `n_bins` is used, then
        `bin_range` is optional. `n_bins` overrides `bin_width`.

        If no options are given, the default is to use 20 bins and the
        range generated by np.histogram.
        """
        # this is to compare whether another histogram had the same setup,
        # and is useful for other programs that want to cache a histogram
        self._inputs = [n_bins, bin_width, bin_range]

        # regularize options
        self.bin_width = bin_width
        self.bin_range = bin_range
        if bin_range is not None:
            max_bin = max(bin_range)
            min_bin = min(bin_range)
            if bin_width is not None:
                self.n_bins = int(math.ceil((max_bin-min_bin)/self.bin_width))
                # if this isn't actually divisible, you'll get one extra bin
            if n_bins is not None:
                self.n_bins = n_bins
                self.bin_width = (max_bin-min_bin)/(self.n_bins)
            self.bins = [min_bin + self.bin_width*i
                         for i in range(self.n_bins+1)]
        else:
            if n_bins is not None:
                self.n_bins = n_bins
                self.bin_width = None
            else:
                self.n_bins = 20  # default
            self.bins = self.n_bins

        try:
            left_bin_edges = (self.bins[0],)
        except TypeError:
            left_bin_edges = None

        super(Histogram, self).__init__(bin_widths=(self.bin_width,),
                                        left_bin_edges=left_bin_edges)

    def empty_copy(self):
        return type(self)(bin_width=self.bin_width, bin_range=self.bin_range)

    def histogram(self, data=None, weights=None):
        """Build the histogram based on `data`.

        Note
        ----
        Calling this with new data overwrites the previous histogram. This
        is the expected behavior; in using this, you should check if the
        histogram parameters have changed from a previous run (using
        `compare_parameters`) and you should be aware whether your data has
        changed. If you want to add data to the histogram, you should use
        `add_data_to_histogram`.
        """
        if self.left_bin_edges is not None:
            return super(Histogram, self).histogram(data, weights)
        if data is not None:
            max_val = max(data)
            min_val = min(data)
            self.bin_width = (max_val-min_val)/self.bins
            self.left_bin_edges = np.array((min_val,))
            self.bin_widths = np.array((self.bin_width,))
        return super(Histogram, self).histogram(data, weights)

    def xvals(self, bin_edge_type="l"):
        int_bins = np.array(list(self._histogram.keys()))[:, 0]
        # always include left_edge_bin as 0 point; always include 0 and
        # greater bin values (but allow negative)
        min_bin = min(min(int_bins), 0)
        n_bins = max(int_bins) - min_bin + 1
        width = self.bin_widths[0]
        left_bins = (self.left_bin_edges[0] + np.arange(n_bins) * width)
        return self._left_edge_to_bin_edge_type(left_bins, width,
                                                bin_edge_type)

    def __call__(self, bin_edge="m"):
        """Return copy of histogram if it has already been built"""
        vals = self.xvals(bin_edge)
        hist = self.histogram()
        bins = sorted(hist.keys())
        min_bin = min(bins[0][0], 0)
        max_bin = bins[-1][0]
        bin_range = range(int(min_bin), int(max_bin)+1)
        hist_list = [hist[(b,)] for b in bin_range]
        return LookupFunction(vals, hist_list)

    def compare_parameters(self, other):
        """Return true if `other` has the same bin parameters as `self`.

        Useful for checking whether a histogram needs to be rebuilt.
        """
        if not super(Histogram, self).compare_parameters(other):
            return False
        if type(other.bins) is not int:
            if type(self.bins) is int:
                return False
            for (t, b) in zip(self.bins, other.bins):
                if t != b:
                    return False
        else:
            return self._inputs == other._inputs
        return True

    def _normalization(self):
        """Return normalization constant (integral over this histogram)."""
        hist = self('l')
        bin_edges = self.xvals('l')
        dx = [bin_edges[i+1] - bin_edges[i] for i in range(len(bin_edges)-1)]
        dx += [dx[-1]]  # assume the "forever" bin is same as last limited
        norm = np.dot(hist.values(), dx)
        return norm

    # Yes, the following could be cached. No, I don't think it is worth it.
    # Keep in mind that we need a separate cache for each one that we build,
    # and that typically it will take almost no time to build one of these
    # (runtime in linear in number of histogram bins). Adding caching
    # complicates the code for no real benefit (you're more likely to suffer
    # from L2 cache misses than to get a speedup).

    def normalized(self, raw_probability=False, bin_edge="m"):
        """Return normalized version of histogram.

        By default (`raw_probability` false), this returns the histogram
        normalized by its integral (according to rectangle-rule
        integration). If `raw_probability` is true, this returns the
        histogram normalized by the sum of the bin counts, with no
        consideration of the bin widths.
        """
        normed_hist = self()  # returns a copy
        nnorm = self._normalization() if not raw_probability else self.count
        norm = 1.0/nnorm
        normed_hist_list = [normed_hist(k) * norm for k in normed_hist.keys()]
        xvals = self.xvals(bin_edge)
        return LookupFunction(xvals, normed_hist_list)

    def cumulative(self, maximum=1.0, bin_edge="r"):
        """Cumulative from the left: number of values less than bin value.

        Use `maximum=None` to get the raw counts.
        """
        cumul_hist = []
        total = 0.0
        hist = self(bin_edge)
        for k in sorted(hist.keys()):
            total += hist(k)
            cumul_hist.append(total)

        cumul_hist = np.array(cumul_hist)
        if total == 0:
            warnings.warn("No non-zero data in the histogram")
        elif maximum is not None:
            cumul_hist *= maximum / total

        xvals = self.xvals(bin_edge)
        return LookupFunction(xvals, cumul_hist)

    def reverse_cumulative(self, maximum=1.0, bin_edge="l"):
        """Cumulative from the right: number of values greater than bin value.

        Use `maximum=None` to get the raw counts.
        """
        cumul_hist = []
        total = 0.0
        hist = self(bin_edge)
        for k in reversed(sorted(hist.keys())):
            total += hist(k)
            cumul_hist.insert(0, total)

        cumul_hist = np.array(cumul_hist)
        if total == 0:
            warnings.warn("No non-zero data in the histogram")
        elif maximum is not None:
            cumul_hist *= maximum / total

        xvals = self.xvals(bin_edge)
        return LookupFunction(xvals, cumul_hist)

    def rebinned(self, scaling):
        """Redistributes histogram bins of width binwidth*scaling

        Exact if scaling is an integer; otherwise uses the assumption that
        original bins were uniformly distributed. Note that the original
        data is not destroyed.
        """
        # TODO
        pass

    def plot_bins(self, scaling=1.0):
        """Bins used in plotting. Scaling useful when plotting `rebinned`"""
        # TODO: add scaling support
        return self.bins[1:]


def histograms_to_pandas_dataframe(hists, fcn="histogram", fcn_args={}):
    """Converts histograms in hists to a pandas data frame"""
    keys = None
    frames = []
    for hist in hists:
        # check that the keys match
        if keys is None:
            keys = hist.xvals()
        for (t, b) in zip(keys, hist.xvals()):
            if t != b:
                raise Warning("Bins don't match up")
        if hist.name is None:
            hist.name = int(hists.index(hist))

        hist_data = {
            "histogram": hist,
            "normalized": hist.normalized,
            "reverse_cumulative": hist.reverse_cumulative,
            "cumulative": hist.cumulative,
            "rebinned": hist.rebinned
        }[fcn](**fcn_args).values()

        bin_edge = {
            "histogram": "m",
            "normalized": "m",
            "reverse_cumulative": "l",
            "cumulative": "r"
        }[fcn]
        xvals = hist.xvals(bin_edge)
        frames.append(pd.DataFrame({hist.name: hist_data}, index=xvals))
    all_frames = pd.concat(frames, axis=1)
    return all_frames.fillna(0.0)


def write_histograms(fname, hists):
    """Writes all histograms in list `hists` to file named `fname`

    If the filename is the empty string, then output is to stdout.
    Assumes that all files should have the same bins.
    """
    pass

# TODO: might as well add a main function to this; read data / weight from
# stdin and output an appropriate histogram depending on some options. Then
# it is both a useful script and a library class!


class Histogrammer(object):
    """
    Basically a dictionary to track what each histogram should be making.
    """
    def __init__(self, f, f_args=None, hist_args=None):
        self.f = f
        self.f_args = f_args
        self._hist_args = hist_args
        self.empty_hist = Histogram(**self._hist_args)

    @property
    def hist_args(self):
        return self._hist_args

    @hist_args.setter
    def hist_args(self, val):
        self._hist_args = val
        self.empty_hist = Histogram(**self._hist_args)


class HistogramPlotter2D(object):
    """
    Convenience tool for plotting 2D histograms and plotting data atop them.

    The difficulty is that matplotlib uses the row/column *numbers* of a
    pandas.DataFrame as the actual internal axis. This class carries all the
    information to properly plot things (even mapping to CVs, if the
    histogram supports that).

    The descriptions below will discuss "real space," "bin space," and
    "frame space." Real space refers to the actual values of the input data.
    Bin space refers to the bins that come out of that for histogramming
    (made into continuous parameters). Frame space is bin space shifted such
    that the lowest bin values are 0.

    Parameters
    ----------
    histogram : :class:`.SparseHistogram`
        input histogram to plot
    normed : bool
        whether to normalize the histogram (using raw_probability=True)
    xticklabels : list of float
        the desired locations for plot xticks, in real space
    yticklabels : list of float
        the desired locations for plot yticks, in real space
    xlim : 2-tuple of (float, float)
        horizontal (x-value) range of (minimum, maximum) bounds for
        displaying the plot
    ylim : 2-tuple of (float, float)
        vertical (y-value) range of (minimum, maximum) bounds for
        displaying the plot
    label_format : string
        Python format-style string for formatting tick labels. Default is
        '{:}'.
    """
    def __init__(self, histogram, normed=True, xticklabels=None,
                 yticklabels=None, xlim=None, ylim=None,
                 label_format="{:}"):
        self.histogram = histogram
        self.normed = normed
        self.xticklabels = xticklabels
        self.yticklabels = yticklabels
        self.xlim = xlim
        self.ylim = ylim
        self.label_format = label_format

        self.xticks_, self.xlim_, self.yticks_, self.ylim_ = self.axes_setup(
            xticklabels, yticklabels, xlim, ylim
        )

    def to_bins(self, alist, dof):
        """Convert real-space values to bin-space values for a given dof

        Parameters
        ----------
        alist : list of float
            input in real-space
        dof : integer (0 or 1)
            degree of freedom; 0 is x, 1 is y

        Returns
        -------
        list of float :
            the outputs in bin-space
        """
        left_edge = self.histogram.left_bin_edges[dof]
        bin_width = self.histogram.bin_widths[dof]
        result = None
        if alist is not None:
            result = (np.asarray(alist) - left_edge) / bin_width
        return result

    def axis_input(self, hist, ticklabels, lims, dof):
        """Get ticks, range, and limits for a given DOF

        Parameters
        ----------
        hist : list of float
            input data from the histogram (bin-space)
        ticklabels : list of float or None
            user-set tick labels for this DOF (real-space)
        lims : 2-tuple (float, float) or None
            user-set plot limits for this DOF
        dof : integer (0 or 1)
            degree of freedom; 0 is x, 1 is y

        Returns
        -------
        ticks_ : list of float or None
            user-set ticks in bin-space
        range_ : list of float
            range for the pandas.DataFrame (bin-space)
        lims_ : 2-tuple (float, float)
            range for plot visualization (bin-space)
        """
        ticks_ = self.to_bins(ticklabels, dof)
        lims_ = self.to_bins(lims, dof)
        ticks = [] if ticks_ is None else list(ticks_)
        lims = [] if lims_ is None else list(lims_)
        range_ = (int(min(list(hist) + ticks + lims)),
                  int(max(list(hist) + ticks + lims)))
        if lims_ is None:
            lims_ = (0, range_[1] - range_[0])
        else:
            lims_ = (lims_[0] - range_[0], lims_[1] - range_[0])
        return (ticks_, range_, lims_)

    def axes_setup(self, xticklabels, yticklabels, xlim, ylim):
        r"""Set up both x-axis and y-axis for plotting.

        Also sets self.xrange\_ and self.yrange\_, which are the (bin-space)
        bounds for the pandas.DataFrame.

        Parameters
        ----------
        xticklabels : list of float
            the desired locations for plot xticks, in real space
        yticklabels : list of float
            the desired locations for plot yticks, in real space
        xlim : 2-tuple of (float, float)
            horizontal (x-value) range of (minimum, maximum) bounds for
            displaying the plot
        ylim : 2-tuple of (float, float)
            vertical (y-value) range of (minimum, maximum) bounds for
            displaying the plot

        Returns
        -------
        xticks_ : list of float or None
            user-set xticks in bin-space
        xlim_ : 2-tuple (float, float)
            range in x for plot visualization (bin-space)
        yticks_ : list of float or None
            user-set yticks in bin-space
        ylim_ : 2-tuple (float, float)
            range in y for plot visualization (bin-space)
        """
        if xticklabels is None:
            xticklabels = self.xticklabels
        if yticklabels is None:
            yticklabels = self.yticklabels
        if xlim is None:
            xlim = self.xlim
        if ylim is None:
            ylim = self.ylim
        x, y = list(zip(*self.histogram._histogram.keys()))
        xticks_, xrange_, xlim_ = self.axis_input(x, xticklabels, xlim, dof=0)
        yticks_, yrange_, ylim_ = self.axis_input(y, yticklabels, ylim, dof=1)
        self.xrange_ = xrange_
        self.yrange_ = yrange_
        return (xticks_, xlim_, yticks_, ylim_)

    def ticks_and_labels(self, ticks, ax, dof):
        """Obtain the plot ticks and tick labels for given dof.

        Parameters
        ----------
        ticks : list of float or None
            user-set input (bin-space) for tick locations
        ax : matplotlib.Axes
            axes from the plot
        dof : integer (0 or 1)
            degree of freedom; 0 is x, 1 is y

        Returns
        -------
        ticks : list of float
            tick locations (bin-space, suitable for matplotlib)
        labels : list of string
            labels for the ticks
        """
        if dof == 0:
            ax_ticks = ax.get_xticks()
            minval = self.xrange_[0]
            bw = self.histogram.bin_widths[0]
            edge = self.histogram.left_bin_edges[0]
        elif dof == 1:
            ax_ticks = ax.get_yticks()
            minval = self.yrange_[0]
            bw = self.histogram.bin_widths[1]
            edge = self.histogram.left_bin_edges[1]
        else:  # pragma: no cover
            raise RuntimeError("Bad DOF: " + str(dof))
        to_val = lambda n: (n + minval) * bw + edge
        ticks = ticks if ticks is not None else ax_ticks
        labels = [self.label_format.format(to_val(n)) for n in ticks]
        return (ticks, labels)

    def plot(self, normed=None, xticklabels=None, yticklabels=None,
             xlim=None, ylim=None, **kwargs):
        """Plot the histogram.

        Parameters
        ----------
        normed : bool
            whether to normalize the histogram (using raw_probability=True)
        xticklabels : list of float
            the desired locations for plot xticks, in real space
        yticklabels : list of float
            the desired locations for plot yticks, in real space
        xlim : 2-tuple of (float, float)
            horizontal (x-value) range of (minimum, maximum) bounds for
            displaying the plot
        ylim : 2-tuple of (float, float)
            vertical (y-value) range of (minimum, maximum) bounds for
            displaying the plot
        kwargs :
            additional arguments to pass to plt.pcolormesh

        Returns
        -------
        PolyCollection :
            return value of plt.pcolormesh
        """
        if normed is None:
            normed = self.normed

        xticks_, xlim_, yticks_, ylim_ = self.axes_setup(
            xticklabels, yticklabels, xlim, ylim
        )

        if normed:
            hist_fcn = self.histogram.normalized(raw_probability=True)
        else:
            hist_fcn = self.histogram()
        df = hist_fcn.df_2d(x_range=self.xrange_, y_range=self.yrange_)
        self.df = df

        mesh = plt.pcolormesh(df.fillna(0.0).transpose(), **kwargs)

        (xticks, xlabels) = self.ticks_and_labels(xticks_, mesh.axes, dof=0)
        (yticks, ylabels) = self.ticks_and_labels(yticks_, mesh.axes, dof=1)

        mesh.axes.set_xticks(xticks)
        mesh.axes.set_yticks(yticks)
        mesh.axes.set_xticklabels(xlabels)
        mesh.axes.set_yticklabels(ylabels)
        plt.xlim(xlim_[0], xlim_[1])
        plt.ylim(ylim_[0], ylim_[1])
        plt.colorbar()
        return mesh

    def plot_trajectory(self, trajectory, *args, **kwargs):
        """Plot a trajectory (or CV trajectory) on the axes.

        Additional arguments pass to plt.plot.

        Parameters
        ----------
        trajectory : :class:`.Trajectory` or list of 2-tuple
            list to plot; paths.Trajectory allowed if the histogram can
            convert it to CVs.
        """
        x, y = list(zip(*self.histogram.map_to_float_bins(trajectory)))
        px = np.asarray(x) - self.xrange_[0]
        py = np.asarray(y) - self.yrange_[0]
        plt.plot(px, py, *args, **kwargs)
