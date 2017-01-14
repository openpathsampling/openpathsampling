import logging
import time

from openpathsampling.pathsimulator import PathSampling, MCStep
import openpathsampling as paths
from openpathsampling.ops_logging import initialization_logging

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


class StepSampling(PathSampling):
    """
    General path sampling code.

    Takes a single move_scheme and generates samples from that, keeping one
    per replica after each move.
    """

    calc_name = "StepSampling"

    def __init__(
            self,
            storage=None,
            steplist=None,
            move_scheme=None,
            sample_set=None,
            initialize=True
    ):
        self.steplist = steplist
        super(StepSampling, self).__init__(storage, move_scheme, sample_set, initialize)

    def save_current_step(self):
        """
        Save the current step to the storage

        """
        if self._current_step is not None:
            if self.storage is not None:
                self.storage.steps.save(self._current_step)
            if self.steplist is not None:
                self.steplist.append(self._current_step)