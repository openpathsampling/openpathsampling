import openpathsampling as paths
from openpathsampling.analysis import SparseHistogram

from collections import Counter
import numpy as np

class PathHistogram(SparseHistogram):
    def __init__(self, left_bin_edges, bin_widths, interpolate=True,
                 per_traj=True):
        super(PathHistogram, self).__init__(left_bin_edges=left_bin_edges, 
                                            bin_widths=bin_widths)
        if interpolate is True:
            interpolate = "subdivide"
        self.interpolate = interpolate
        self.per_traj = per_traj

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

        # otherwise, use one of the interpolation algos to find bins
        if self.interpolate == "subdivide":
            bin_list = self.subdivide_interpolation(start_pt=old_pt,
                                                    end_pt=new_pt,
                                                    start_bin=old_bin,
                                                    end_bin=new_bin)
        return list(set(bin_list) - set([old_bin]))

    def subdivide_interpolation(self, start_pt, end_pt, start_bin, end_bin):
        mid_pt = start_pt + 0.5*(np.asarray(end_pt) - np.asarray(start_pt))
        mid_bin = self.map_to_bins(mid_pt)
        if mid_bin == start_bin:
            return self.subdivide_interpolation(start_pt=mid_pt,
                                                end_pt=end_pt,
                                                start_bin=mid_bin,
                                                end_bin=end_bin)
        elif mid_bin == end_bin:
            return self.subdivide_interpolation(start_pt=start_pt,
                                                end_pt=mid_pt,
                                                start_bin=start_bin,
                                                end_bin=mid_bin)

        manhattan_dist_start = sum(abs(np.asarray(mid_bin) - 
                                       np.asarray(start_bin)))
        manhattan_dist_end = sum(abs(np.asarray(end_bin) - 
                                     np.asarray(mid_bin)))
        if manhattan_dist_start == 1 and manhattan_dist_end == 1:
            return [start_bin, mid_bin, end_bin]
        elif manhattan_dist_start == 1:
            return ([start_bin] + 
                    self.subdivide_interpolation(start_pt=mid_pt,
                                                 end_pt=end_pt,
                                                 start_bin=mid_bin,
                                                 end_bin=end_bin))
        elif manhattan_dist_end == 1:
            return ([end_bin] +
                    self.subdivide_interpolation(start_pt=start_pt,
                                                 end_pt=mid_pt,
                                                 start_bin=start_bin,
                                                 end_bin=mid_bin))
        else:
            start_side = self.subdivide_interpolation(start_pt=start_pt,
                                                      end_pt=mid_pt,
                                                      start_bin=start_bin,
                                                      end_bin=mid_bin)
            end_side = self.subdivide_interpolation(start_pt=mid_pt,
                                                    end_pt=end_pt,
                                                    start_bin=mid_bin,
                                                    end_bin=end_bin)
            return start_side + end_side


    def add_trajectory(self, traj):
        # make a list of every bin visited, possibly interpolating gaps
        bin_list = [self.map_to_bins(traj[0])]
        for fnum in range(len(traj)-1):
            if self.interpolate:
                bin_list += self.interpolated_bins(traj[fnum], traj[fnum+1])
            else:
                bin_list += [self.map_to_bins(traj[fnum+1])]

        local_hist = Counter(bin_list)
        if self.per_traj:
            # keys only exist once, so the counter gives 1 if key present
            local_hist = Counter(local_hist.keys())
        if self._histogram is None:
            self._histogram = Counter({})
        self._histogram += local_hist


class PathDensityHistogram(SparseHistogram):
    def __init__(self, cvs, left_bin_edges, bin_widths):
        pass
