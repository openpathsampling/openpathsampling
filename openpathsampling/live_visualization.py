import matplotlib
import matplotlib.pyplot as plt

import openpathsampling as paths

class LiveVisualization(object):
    def __init__(self, network, cv_x, cv_y, xlim, ylim, output_directory=None):
        self.network = network
        self.cv_x = cv_x
        self.cv_y = cv_y
        self.xlim = xlim
        self.ylim = ylim
        self.output_directory = output_directory
        self.background = None
        self._save_bg_axes = None

        self.extent = self.xlim + self.ylim

        self.add_lines = [] # keeps state

        default_colors = {0 : 'r', 1 : 'b', 2 : 'g', 3 : 'y', 4 : 'm', 5 : 'c'}
        self.ensemble_colors = {}
        for transition in self.network.sampling_transitions:
            for i in range(len(transition.ensembles)):
                ens = transition.ensembles[i]
                self.ensemble_colors[ens] = default_colors[i % 6]
        for special_type in self.network.special_ensembles:
            for ens in self.network.special_ensembles[special_type].keys():
                self.ensemble_colors[ens] = 'k'
            

    def draw(self, mcstep):
        self.fig, self.ax = plt.subplots()
        if self.background is not None:
            if self._save_bg_axes is None:
                self._save_bg_axes = self.background.axes
            self.fig = self.background
            for ax in self.fig.axes:
                self.fig.delaxes(ax)
            for ax in self._save_bg_axes:
                self.fig.add_axes(ax)
            self.ax = self.fig.axes[0].twinx()
            self.ax.cla()
            #ax = self.fig.axes[0]
        #else:
            #self.fig, ax = plt.subplots()
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)
        # draw the background
        # draw the active trajectories
        for sample in mcstep.active:
            #lines.append(matplotlib.lines.lines2D(
            self.ax.plot(
                self.cv_x(sample.trajectory), 
                self.cv_y(sample.trajectory), 
                linewidth=2, zorder=2,
                color=self.ensemble_colors[sample.ensemble]
            )#)
            # draw arrowheads at the end of each active

        # draw the trials
        for trial in mcstep.change.trials:
            #lines.append(matplotlib.lines.lines2D(
            self.ax.plot(
                self.cv_x(trial.trajectory), 
                self.cv_y(trial.trajectory), 
                linewidth=1.0, zorder=3,
                color=self.ensemble_colors[trial.ensemble]
            )#)
            # draw arrowheads at the end of each trial
        # decorate by trial type

        return self.fig


    def draw_ipynb(self, mcstep):
        try:
            import IPython.display
        except ImportError:
            pass
        else:
            IPython.display.clear_output(wait=True)
            fig = self.draw(mcstep)
            IPython.display.display(fig);
            plt.close() # prevents crap in the output

    def draw_png(self, mcstep):
        pass
