from functools import partial
from .options import canonicalize_mover
from .plotter import PathTreePlotter

# In addition to the standard `draw_trajectory` and `draw_connector` methods
# required of all PathTree plot styles, plot styles based on matplotlib
# should have the attribute `ax`, which returns the `Axes` object they plot
# into.



class MPLPlottingStyle(PathTreePlotter):
    """Abstract base class for matplotlib-based path tree plotting styles.

    Parameters
    ----------
    ax : :class:`.matplotlib.Axes`
        Axes object to plot into
    """
    def __init__(self, ax=None, options=None):
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise RuntimeError("matplotlib is not installed")

        super().__init__(options)
        if ax is None:
            _, ax = plt.subplots()
        self.ax = ax

    # @staticmethod
    # def step_basics(step, options):
    #     mover = canonicalize_mover(step.mover)
    #     mover_options = options.movers[mover]
    #     color = mover_options.color
    #     plot_segments = mover_options.get_left_right(step)
    #     return mover, color, plot_segments

    def draw_trajectory(self, row, step, options):
        raise NotImplementedError()

    def draw_connector(self, x, bottom, top, step):
        mover = canonicalize_mover(step.mover)
        if mover in self.CONNECTOR_FIXES:
            x = self.CONNECTOR_FIXES[mover](x)

        self.ax.plot([x,x], [-bottom, -top], color='k')


### ThinLines plotting style ###############################################

class ThinLines(MPLPlottingStyle):
    """ThinLines plotting style for matplotlib-based path trees.

    This is a plotting style based on using `ax.plot` and with the style
    choice that the shooting point is always represented at the middle of
    its frame. This needs to us filling half of the shooting point frame so
    that the trajectory line intersects the shooting point connector. This
    is accomplished by the ``thin_lines_fix`` function.

    """
    PLOT_TYPE = "trajectory"
    def draw_trajectory(self, row, step):
        mover = canonicalize_mover(step.mover)
        plot_segments, color = self.get_step_plot_details(step)

        for left, right in plot_segments:
            self.ax.plot([left - 0.5, right - 0.5], [-row] * 2, lw=3,
                         color=color)



### Blocks plotting style ##################################################

from matplotlib.patches import Rectangle

class Blocks(MPLPlottingStyle):
    """Blocks plotting style for matplotlib-based path trees.

    The Blocks plotting style uses rectangular patches for each snapshot.
    """
    PLOT_TYPE = "snapshot"
    def __init__(self, ax=None, options=None, *, block_size=(0.8, 0.6)):
        super().__init__(ax, options)
        self.block_size = block_size

    def draw_trajectory(self, row, step):
        """
        Parameters
        ----------
        row : int
            the row number (counting from the top) for this step
        """
        mover = canonicalize_mover(step.mover)
        plot_segments, color = self.get_step_plot_details(step)
        # first we draw the background line
        for left, right in plot_segments:
            self.ax.plot([left - 0.5, right - 0.5], [-row]*2, color=color)

        # now we draw the boxes
        for left, right in plot_segments:
            width, height = self.block_size
            y = - row - height / 2
            for snap in range(left, right):
                self.ax.add_patch(
                    Rectangle((snap - width / 2, y), width, height,
                              color=color)
                )
