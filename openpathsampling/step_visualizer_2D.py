import matplotlib
import matplotlib.pyplot as plt

import openpathsampling as paths
import logging
logger = logging.getLogger(__name__)

class StepVisualizer2D(object):
    def __init__(self, network, cv_x, cv_y, xlim, ylim, output_directory=None):
        self.network = network
        self.cv_x = cv_x
        self.cv_y = cv_y
        self.xlim = xlim
        self.ylim = ylim
        self.output_directory = output_directory
        is_interactive = matplotlib.is_interactive()
        self.background = None
        self._save_bg_axes = None

        self.fig = None

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

    def draw_arrowhead(self, sample, accepted=True):
        if accepted:
            face_color = self.ensemble_colors[sample.ensemble]
        else:
            face_color = 'none'
        # for now, our "arrowheads" are just circles
        frame = sample.trajectory[-1]
        self.ax.plot(
            self.cv_x(frame), self.cv_y(frame), marker='o',
            markeredgecolor=self.ensemble_colors[sample.ensemble],
            markerfacecolor=face_color,
            zorder=5
        )

    def draw_samples(self, samples, accepted=True):
        self.draw_background()
        for sample in samples:
            if sample.ensemble in self.ensemble_colors:
                self.ax.plot(
                    self.cv_x(sample.trajectory),
                    self.cv_y(sample.trajectory),
                    linewidth=2, zorder=2,
                    color=self.ensemble_colors[sample.ensemble]
                )
                # draw arrowheads at the end of each active
                self.draw_arrowhead(sample, accepted=accepted)
        return self.fig

    def draw_trials(self, change):
        if self.fig != self.background:
            self.draw_background()
        for trial in change.trials:
            self.ax.plot(
                self.cv_x(trial.trajectory),
                self.cv_y(trial.trajectory),
                linewidth=1.0, zorder=3,
                color=self.ensemble_colors[trial.ensemble]
            )
            # draw arrowheads at the end of each trial
            self.draw_arrowhead(trial, accepted=False)
        return self.fig

    def draw_background(self):
        # draw the background
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
        else:
            # create empty figure and axis
            self.fig, self.ax = plt.subplots()
            # set the empty axis and figure as our new 'background'
            # this avoids flickering if we always recreate the figure
            self.background = self.fig
            self._save_bg_axes = self.background.axes
            # create a copy axis to draw on
            self.ax = self.fig.axes[0].twinx()

        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)
        return

    def draw(self, mcstep):
        self.draw_samples(mcstep.active)
        self.draw_trials(mcstep.change)
        return self.fig

    def draw_ipynb(self, mcstep):
        try:
            import IPython.display
        except ImportError:
            logger.info("Not in IPython")
        else:
            IPython.display.clear_output(wait=True)
            fig = self.draw(mcstep)
            IPython.display.display(fig);
            plt.close() # prevents crap in the output

    def draw_png(self, mcstep):
        pass
