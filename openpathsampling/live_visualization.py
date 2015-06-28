import matplotlib
import matplotlib.pyplot as plt

import openpathsampling as paths

class LiveVisualization(object):
    def __init__(self, cv_x, cv_y, range_x, range_y, output_directory=None):
        self.cv_x = cv_x
        self.cv_y = cv_y
        self.range_x = range_x
        self.range_y = range_y
        self.output_directory = output_directory
        self.extent = [range_x[0], range_x[1], range_y[0], range_y[1]]

    def draw(self, mcstep):
        try:
            import IPython.display
        except ImportError:
            pass
        else:
            IPython.display.clear_output(wait=True)
            fig, ax = plt.subplots()
            # draw the trajectories
            for traj in mcstep.active:
                ax.plot(self.cv_x(traj), self.cv_y(traj));

            IPython.display.display(fig);
            plt.close() # prevents crap in the output
