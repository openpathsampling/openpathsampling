import openpathsampling as paths
import math

import numpy as np

import matplotlib.pyplot as plt


class CVSphere(paths.Volume):
    """
    Defines a sphere in multi-CV space with center and distance
    """

    def __init__(self, cvs, center, radius):
        super(CVSphere, self).__init__()
        self.cvs = cvs
        self.center = center
        self.radius = radius

        assert (len(cvs) == len(center) == len(radius))

    def __call__(self, snapshot):
        return math.sqrt(sum(
            map(
                lambda cv: cv(snapshot) ** 2
            ), self.cvs
        ))

    def __and__(self, other):
        if isinstance(other, paths.EmptyVolume):
            return self
        elif isinstance(other, paths.FullVolume):
            return other
        elif isinstance(other, CVSphere):
            dc = np.linalg.norm(np.array(self.center) - np.array(other.center))

            # use triangle inequality
            if self.radius >= dc + other.radius:
                # other is completely in self
                return self
            elif other.radius >= dc + self.radius:
                # self is completely in other
                return other

        return paths.UnionVolume(
            self, other
        )


class TwoCVSpherePlot(object):
    def __init__(
            self, cvs, states, state_centers,
            interface_levels, ranges=None):
        self.cvs = cvs
        self.states = states
        self.state_centers = state_centers
        self.interface_levels = interface_levels
        self._ax1 = 0
        self._ax2 = 1
        self.figsize = (6, 6)
        self.periodic = [math.pi] * len(cvs)
        self.zoom = 180 / math.pi
        if ranges is None:
            self.ranges = ((-180, 180), (-180, 180))
        else:
            self.ranges = ranges

        self.markersize = 2
        self.color_fnc = lambda x: (x, x, 0.6)
        self.color_fnc = lambda x: (x * 0.5 + 0.4, 0.5 * x + 0.4, 1 * x, 1.0)

    def select_axis(self, ax1, ax2):
        self._ax1 = ax1
        self._ax2 = ax2

    def new(self, figsize=None):
        if figsize is None:
            figsize = self.figsize

        plt.figure(figsize=figsize)

    def main(self):
        n_states = len(self.states)

        centers = self.state_centers
        levels = self.interface_levels
        labels = [state.name[0] for state in self.states]
        periodic = (self.periodic[self._ax1], self.periodic[self._ax2])

        mirror = [
            [-1, 0, 1] if p is not None else [0]
            for p in periodic
            ]

        # replace None with zero
        periodic = [p or 0 for p in periodic]

        plt.plot(
            [x[self._ax1] for x in centers],
            [x[self._ax2] for x in centers],
            'ko',
            markersize=self.markersize)

        fig = plt.gcf()

        all_levels = sorted(
            list(set(
                sum(levels, [])
            )),
            reverse=True
        ) + [0]

        plt.xlabel(self.cvs[self._ax1].name)
        plt.ylabel(self.cvs[self._ax2].name)

        max_level = max(all_levels)
        zoom = self.zoom

        for level in all_levels:
            for colored in [True, False]:
                for state in range(n_states):
                    center = centers[state]
                    center = (center[self._ax1], center[self._ax2])
                    name = labels[state]

                    if level == 0:
                        plt.annotate(
                            name,
                            xy=center,
                            xytext=(center[0] + 10 + 1, center[1] - 1),
                            fontsize=20,
                            color='k'
                        )
                        plt.annotate(
                            name,
                            xy=center,
                            xytext=(center[0] + 10, center[1]),
                            fontsize=20,
                            color='w'
                        )

                    if level in levels[state]:
                        for xp in mirror[0]:
                            for yp in mirror[1]:
                                if colored:
                                    circle = plt.Circle(
                                        (center[0] + xp * periodic[
                                            0] * zoom * 2,
                                         center[1] + yp * periodic[
                                             1] * zoom * 2),
                                        level,
                                        color='w'
                                    )
                                    fig.gca().add_artist(circle)
                                else:
                                    l = 1.0 * level / max_level
                                    circle = plt.Circle(
                                        (center[0] + xp * periodic[
                                            0] * zoom * 2,
                                         center[1] + yp * periodic[
                                             1] * zoom * 2),
                                        level - 1,
                                        color=self.color_fnc(l)
                                    )
                                    fig.gca().add_artist(circle)

        # plt.axis((-180,180,-180,180))

        plt.axis('equal')
        plt.xlim(*self.ranges[0])
        plt.ylim(*self.ranges[1])

    def _cvlines(self, snapshots):
        cvs = self.cvs
        all_points = [cv(snapshots) for cv in cvs]
        ret = []
        first = 0
        if len(snapshots) > 1:
            for d in range(1, len(snapshots)):
                flip = False
                for c in range(len(cvs)):
                    if self.periodic[c] is not None and self._periodicflip(
                            all_points[c][d],
                            all_points[c][d - 1],
                            self.periodic[c]
                    ):
                        flip = True

                if flip:
                    ret.append(
                        [all_points[c][first:d] for c in range(len(cvs))])
                    first = d

            ret.append([all_points[c][first:d + 1] for c in range(len(cvs))])

        return ret

    @staticmethod
    def _periodicflip(val1, val2, period):
        return (period ** 2 - (val1 - val2) ** 2) < (val1 - val2) ** 2

    def add_trajectory(self, trajectory, line=True, points=True):
        angles = self._cvlines(trajectory)
        zoom = self.zoom

        for angle in angles:
            if points:
                plt.plot(
                    zoom * np.array(angle[self._ax1])[:],
                    zoom * np.array(angle[self._ax2])[:],
                    'ko',
                    markersize=self.markersize,
                    linewidth=0.5)
            if line:
                plt.plot(
                    zoom * np.array(angle[self._ax1])[:],
                    zoom * np.array(angle[self._ax2])[:],
                    'k-',
                    markersize=self.markersize,
                    linewidth=0.5)

    def add_snapshot(self, snapshot, label=None):
        zoom = self.zoom
        angle = [cv(snapshot) for cv in self.cvs]

        x = zoom * np.array(angle[self._ax1])
        y = zoom * np.array(angle[self._ax2])

        plt.plot(
            x, y,
            'w+',
            mew=5, ms=14)

        plt.plot(
            x, y,
            'k+',
            mew=3, ms=12)

        if label is not None:
            plt.annotate(
                label,
                xy=(x, y),
                xytext=(x + 6, y + 4),
                fontsize=12,
                color='w'
            )
            plt.annotate(
                label,
                xy=(x, y),
                xytext=(x + 5, y + 5),
                fontsize=12,
                color='k'
            )