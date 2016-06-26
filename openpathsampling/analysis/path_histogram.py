import openpathsampling as paths
from openpathsampling.analysis import SparseHistogram

from collections import Counter
import numpy as np

class PathHistogram(SparseHistogram):
    def __init__(self, left_bin_edges, bin_widths):
        super(PathHistogram, self).__init__(left_bin_edges=left_bin_edges, 
                                            bin_widths=bin_widths)

    def interpolated_bins(self, old_pt, new_pt):
        # bins for this do not include the bin of the old point
        old_bin = self.map_to_bins(old_pt)
        new_bin = self.map_to_bins(new_pt)
        abs_dx = abs(np.asarray(new_bin) - np.asarray(old_bin))
        manhattan_distance = sum(abs_dx)
        bin_list = [new_bin]
        # if the manhattan_distance is 1, we're adjacent
        if manhattan_distance == 1:
            return bin_list
        # if the abs_dx is the same for each, we should check if we're
        # unluckily on the diagonal
        if np.all(abs_dx == abs_dx[0]):
            pass # test that we aren't on diag

        # otherwise, subdivide until adjacent
        return bin_list

    def add_trajectory(self, traj, interpolate=True, per_traj=True):
        # make a list of every bin visited, possibly interpolating gaps
        bin_list = [self.map_to_bins(traj[0])]
        for fnum in range(len(traj)-1):
            if interpolate:
                bin_list += self.interpolated_bins(traj[fnum], traj[fnum+1])
            else:
                bin_list += [self.map_to_bins(traj[fnum+1])]

        local_hist = Counter(bin_list)
        if per_traj:
            # keys only exist once, so the counter gives 1 if key present
            local_hist = Counter(local_hist.keys())
        if self._histogram is None:
            self._histogram = Counter({})
        self._histogram += local_hist


class PathDensityHistogram(SparseHistogram):
    def __init__(self, cvs, left_bin_edges, bin_widths):
        pass
