"""
ABCs to clarify types of plotting styles, to make it easier to implement
others in the future.
"""

from .options import canonicalize_mover, create_default_options

class PathTreePlotter(object):
    """
    """
    @property
    def PLOT_TYPE(self):
        raise NotImplementedError("Subclasses must be either 'snapshot' or "
                                  "'trajectory' for PLOT_TYPE")
    def __init__(self, options=None):
        if options is None:
            options = create_default_options()

        self.options = options

    def get_step_plot_details(self, step):
        mover = canonicalize_mover(step.mover)
        mover_options = self.options.movers[mover]
        color = mover_options.color
        plot_segments = mover_options.get_left_right(step)
        return plot_segments, color

    def draw_trajectory(self, row, step):
        raise NotImplementedError()

    def draw_connector(self, x, bottom, top, step):
        raise NotImplementedError()

    def draw_trajectories(self, steps, block_offset=0):
        mccycle_mapping = {}
        for row, step in enumerate(self.options.filter_tree_steps(steps)):
            mccycle_mapping[step.mccycle] = row

            self.draw_trajectory(row, step)

            if step.connector is not None:
                x = step.offset + step.connector.frame
                bottom = mccycle_mapping[step.mccycle]
                top = mccycle_mapping[step.connector.mccycle]
                self.draw_connector(x, bottom, top, step)

    def draw(self, path_tree):
        self.draw_trajectories(path_tree.path_tree_steps)
