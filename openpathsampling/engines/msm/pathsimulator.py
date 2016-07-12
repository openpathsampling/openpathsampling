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

    calc_name = "PathSampling"
    def __init__(
            self,
            storage=None,
            steplist=None,
            move_scheme=None,
            globalstate=None
    ):
        """
        Parameters
        ----------
        storage : openpathsampling.storage.Storage
            the storage where all results should be stored in
        move_scheme : openpathsampling.MoveScheme
            the move scheme used for the pathsampling cycle
        globalstate : openpathsampling.SampleSet
            the initial SampleSet for the Simulator
        """
        super(PathSampling, self).__init__(storage)
        self.move_scheme = move_scheme
        self.root_mover = move_scheme.move_decision_tree()

        self.steplist = steplist

        samples = []
        if globalstate is not None:
            for sample in globalstate:
                samples.append(sample.copy_reset())

        self.globalstate = paths.SampleSet(samples)
        self.root = self.globalstate

        initialization_logging(init_log, self,
                               ['move_scheme', 'globalstate'])
        self.live_visualization = None
        self.visualize_frequency = 1
        self._mover = paths.PathSimulatorMover(self.root_mover, self)

    def run(self, nsteps):
        mcstep = None

        cvs = list()
        n_samples = 0

        if self.storage is not None:
            n_samples = len(self.storage.snapshots)
            cvs = list(self.storage.cvs)

        if self.step == 0:
            if self.storage is not None:
                self.storage.save(self.move_scheme)
            mcstep = self.save_initial()

            if self.steplist is not None:
                self.steplist.append(mcstep)

        for nn in range(nsteps):
            self.step += 1
            logger.info("Beginning MC cycle " + str(self.step))
            refresh=True
            if self.step % self.visualize_frequency == 0:
                # do we visualize this step?
                if self.live_visualization is not None and mcstep is not None:
                    # do we visualize at all?
                    self.live_visualization.draw_ipynb(mcstep)
                    refresh = False

                paths.tools.refresh_output(
                    "Working on Monte Carlo cycle number " + str(self.step)
                    + ".\n",
                    refresh=refresh
                )

            time_start = time.time()
            movepath = self._mover.move(self.globalstate, step=self.step)
            samples = movepath.results
            new_sampleset = self.globalstate.apply_samples(samples)
            time_elapsed = time.time() - time_start

            # TODO: we can save this with the MC steps for timing? The bit
            # below works, but is only a temporary hack
            setattr(movepath.details, "timing", time_elapsed)

            mcstep = MCStep(
                simulation=self,
                mccycle=self.step,
                previous=self.globalstate,
                active=new_sampleset,
                change=movepath
            )

            if self.storage is not None:
                for cv in cvs:
                    n_len = len(self.storage.snapshots)
                    cv(self.storage.snapshots[n_samples:n_len])
                    n_samples = n_len

                self.storage.steps.save(mcstep)

            if self.steplist is not None:
                self.steplist.append(mcstep)

            if self.step % self.save_frequency == 0:
                self.globalstate.sanity_check()
                self.sync_storage()

            self.globalstate = new_sampleset

        self.sync_storage()

        if self.live_visualization is not None and mcstep is not None:
            self.live_visualization.draw_ipynb(mcstep)
        paths.tools.refresh_output(
            "DONE! Completed " + str(self.step) + " Monte Carlo cycles.\n",
            refresh=False
        )