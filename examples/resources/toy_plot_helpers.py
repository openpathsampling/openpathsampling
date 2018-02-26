import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np

from openpathsampling.engines.toy import Snapshot


# A little class we use for visualizing these 2D PESs
class CallablePES(object):
    def __init__(self, pes):
        self.pes = pes

    def __call__(self, x, y):
        self.positions = [x, y]
        return self.pes.V(self)

class CallableVolume(object):
    def __init__(self, vol):
        self.vol = vol

    def __call__(self, x, y):
        snapshot = Snapshot(coordinates=np.array([[x,y,0.0]]))
        return 1.0 if self.vol(snapshot) else 0.0

class ToyPlot(object):
    def __init__(self):
        range_x = np.arange(-1.1, 1.1, 0.01)
        range_y = np.arange(-1.1, 1.1, 0.01)
        self.extent = [range_x[0], range_x[-1], range_y[0], range_y[-1]]
        self.X, self.Y = np.meshgrid(range_x, range_y)
        pylab.rcParams['figure.figsize'] = 9, 6
        self.repcolordict = {0 : 'k-', 1 : 'r-', 2 : 'g-', 3 : 'b-',
                             4 : 'r-'}

        self.contour_range = np.arange(0.0, 1.5, 0.1)

        self._states = None
        self._pes = None
        self._interfaces = None
        self._initcond = None


    def add_pes(self, pes):
        if self._pes is None:
            self._pes = np.vectorize(CallablePES(pes))(self.X, self.Y)

    def add_states(self, states):
        if self._states is None:
            state = states[0]
            self._states = np.vectorize(CallableVolume(state))(self.X, -self.Y)
            for state in states[1:]:
                self._states += np.vectorize(CallableVolume(state))(self.X, -self.Y)

    def add_interfaces(self, ifaces):
        if self._interfaces is None:
            self._interfaces = []
            for iface in ifaces:
                self._interfaces.append(
                    np.vectorize(CallableVolume(iface))(self.X,self.Y)
                )

    def add_initial_condition(self, initcond):
        self._initcond = initcond

    def plot_pes_initcond(self, trajectories):
        fig, ax = plt.subplots()
        if self._pes is not None:
            plt.contour(self.X, self.Y, self._pes,
                        levels=np.arange(0.0, 1.5, 0.1), colors='k')
        if self._initcond is not None:
            ax.plot(self._initcond.coordinates[0,0],
                    self._initcond.coordinates[0,1],
                    'ro', zorder=3)
        for traj in trajectories:
            plt.plot(traj.coordinates()[:,0,0], traj.coordinates()[:,0,1],
                     self.repcolordict[trajectories.index(traj)],
                     zorder=2)

    def plot(self, trajectories=[], bold=[]):
        fig, ax = plt.subplots()
        if self._states is not None:
            plt.imshow(self._states, extent=self.extent, cmap="Blues",
                       interpolation='nearest', vmin=0.0, vmax=2.0,
                       aspect='auto')
        if self._pes is not None:
            plt.contour(self.X, self.Y, self._pes,
                        levels=self.contour_range, colors='k')
        if self._interfaces is not None:
            for iface in self._interfaces:
                plt.contour(self.X, self.Y, iface,
                            colors='r', interpolation='none', levels=[0.5])
        if self._initcond is not None:
            ax.plot(self._initcond.coordinates[0,0],
                    self._initcond.coordinates[0,1],
                    'ro', zorder=3)
        for traj in bold:
            plt.plot(traj.xyz[:,0,0], traj.xyz[:,0,1],
                     self.repcolordict[bold.index(traj)], linewidth=2,
                    zorder=1)
        for traj in trajectories:
            plt.plot(traj.xyz[:,0,0], traj.xyz[:,0,1],
                     self.repcolordict[trajectories.index(traj) % 5],
                     zorder=2)
        return fig

    def reset(self):
        self._pes = None
        self._interfaces = None
        self._initcond = None
        self._states = None


