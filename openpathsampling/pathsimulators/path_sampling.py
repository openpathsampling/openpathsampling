import time
import logging
import os

import openpathsampling as paths
from .path_simulator import PathSimulator, MCStep
from ..ops_logging import initialization_logging
from openpathsampling.beta import hooks


logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')

class PathSampling(PathSimulator):
    """
    General path sampling code.

    Takes a single move_scheme and generates samples from that, keeping one
    per replica after each move.
    """

    calc_name = "PathSampling"

    def __init__(self, storage, move_scheme=None, sample_set=None,
                 initialize=True):
        """
        Parameters
        ----------
        storage : :class:`openpathsampling.storage.Storage`
            the storage where all results should be stored in
        move_scheme : :class:`openpathsampling.MoveScheme`
            the move scheme used for the pathsampling cycle
        sample_set : :class:`openpathsampling.SampleSet`
            the initial SampleSet for the Simulator
        initialize : bool
            if `False` the new PathSimulator will continue at the step and
            not create a new SampleSet object to cut the connection to previous
            steps
        """
        super(PathSampling, self).__init__(storage)
        self.move_scheme = move_scheme
        if move_scheme is not None:
            self.root_mover = move_scheme.move_decision_tree()
            self._mover = paths.PathSimulatorMover(self.root_mover, self)
        else:
            self.root_mover = None
            self._mover = None

        initialization_logging(init_log, self,
                               ['move_scheme', 'sample_set'])
        # TODO: move this to PathSamplingOutputHook defaults?
        # i.e. do not attach a hook on default and use status_update_freq=1
        # as default for PathSamplingOutputHook,
        # let the user add/modify hook if needed?
        # this would also 'fix' the fact that you can only set these two
        # once and before starting any simulation
        # because the hooks read them only before_simulation and
        # if their values for them are None
        self.live_visualizer = None
        self.status_update_frequency = 1

        if initialize:
            # NOTE: why aren't we using save_initial_step here?
            samples = []
            if sample_set is not None:
                for sample in sample_set:
                    samples.append(sample.copy_reset())

            self.sample_set = paths.SampleSet(samples)

            mcstep = MCStep(
                simulation=self,
                mccycle=self.step,
                active=self.sample_set,
                change=paths.AcceptedSampleMoveChange(self.sample_set.samples)
            )

            self._current_step = mcstep

        else:
            self.sample_set = sample_set
            self._current_step = None

        self.root = self.sample_set

        if self.storage is not None:
            template_trajectory = self.sample_set.samples[0].trajectory
            self.storage.save(template_trajectory)
            self.storage.save([self.move_scheme, self.root_mover,
                               self._mover])
            self.save_current_step()

    def to_dict(self):
        return {
            'root': self.root,
            'move_scheme': self.move_scheme,
            'root_mover': self.root_mover,
        }

    @classmethod
    def from_dict(cls, dct):

        # create empty object
        obj = cls(None)

        # and correct the content
        obj.move_scheme = dct['move_scheme']
        obj.root = dct['root']
        obj.root_mover = dct['root_mover']

        obj._mover = paths.PathSimulatorMover(obj.root_mover, obj)

        return obj

    def attach_default_hooks(self):
        self.attach_hook(hooks.StorageHook())
        self.attach_hook(hooks.PathSamplingOutputHook())
        # TODO: or should we make self.live_visualizer a property?
        # and only attach the visualizerhook on 'request', i.e. by setting it?
        self.attach_hook(hooks.LiveVisualizerHook())

    @property
    def current_step(self):
        return self._current_step

    def save_current_step(self):
        """
        Save the current step to the storage

        """
        if self.storage is not None and self._current_step is not None:
            try:
                # new storage does a stash here, not a save
                self.storage.stash(self._current_step)
            except AttributeError:
                self.storage.steps.save(self._current_step)

    @classmethod
    def from_step(cls, storage, step, initialize=True):
        """

        Parameters
        ----------
        storage : :class:`openpathsampling.storage.Storage`
            the storage to be used to hold the simulation results
        step : :class:`openpathsampling.MCStep`
            the step used to fill the initial parameters
        initialize : bool
            if `False` the new PathSimulator will continue at the given step and
            not create a new SampleSet object to cut the connection to previous
            steps.

        Returns
        -------
        :class:`openpathsampling.PathSampling`
            the new simulator object
        """
        obj = cls(
            storage,
            step.simulation.move_scheme,
            step.sample_set,
            initialize=initialize
        )

        return obj

    def restart_at_step(self, step, storage=None):
        """
        Continue with a loaded pathsampling at a given step

        Notes
        -----
        You can only continue from a step that is compatible in the sense
        that it was previously generated from the pathsampling instance.

        If you want to switch the move scheme you need to create a new
        pathsampling instance. You can do so with the constructor or using
        the classmethod `from_step` which simplifies the setup process

        Parameters
        ----------
        step : :class:`MCStep`
            the step to be continued from. You are always free to chose any step
            which can be used to fork a simulation but for analysis you may
            only use one path of steps.
        storage : :class:`Storage`
            If given this will change the storage used to store the generated
            steps

        """
        if step.simulation is not self:
            raise RuntimeWarning(
                'Trying to continue from other step. Please use the '
                '`.from_step` method to create a new PathSampling object '
                'instead.')

        if storage is not None:
            self.storage = storage

        self.step = step.mccycle
        self.sample_set = step.active
        self.root = step.simulation.root

        self._current_step = step

    def run_until(self, n_steps):
        # if self.storage is not None:
        #     if len(self.storage.steps) > 0:
        #         self.step = len(self.storage.steps)
        n_steps_to_run = n_steps - self.step
        self.run(n_steps_to_run)

    def run_until_decorrelated(self, time_reversal=True):
        """Run until all trajectories are decorrelated.

        This runs until all the replicas in ``self.sample_set`` have
        decorrelated from their initial conditions. "Decorrelated" here is
        meant in the sense commonly used in one-way shooting: this runs
        until no configurations from the original trajectories remain.
        """
        originals = {s.replica: s.trajectory for s in self.sample_set}
        current = self.sample_set

        # cache the output stream; force the primary `run` method to not
        # output anything
        original_output_stream = self.output_stream
        self.output_stream = open(os.devnull, 'w')

        def n_correlated(sample_set, originals):
            return sum([originals[r].is_correlated(sample_set[r],
                                                   time_reversal)
                        for r in originals])

        original_output_stream.write("Decorrelating trajectories....\n")
        to_decorrelate = n_correlated(self.sample_set, originals)
        # walrus in py38!
        while to_decorrelate:
            out_str = "Step {}: {} of {} trajectories still correlated\n"
            paths.tools.refresh_output(
                out_str.format(self.step + 1, to_decorrelate, len(originals)),
                refresh=False,
                output_stream=original_output_stream
            )
            self.run(1)
            to_decorrelate = n_correlated(self.sample_set, originals)

        paths.tools.refresh_output(
            "Step {}: All trajectories decorrelated!\n".format(self.step+1),
            refresh=False,
            output_stream=original_output_stream
        )

        self.output_stream = original_output_stream

    def run(self, n_steps):
        mcstep = None
        initial_time = time.time()
        hook_state = None
        self.run_hooks('before_simulation', sim=self)
        for nn in range(n_steps):
            # bookkeeping and before_step hooks
            self.step += 1
            logger.info("Beginning MC cycle " + str(self.step))
            step_number = self.step
            time_start = time.time()
            elapsed_total = time.time() - initial_time
            # TODO: anything else we could want in the step info?
            step_info = (nn, n_steps, elapsed_total)
            self.run_hooks('before_step', sim=self, step_number=step_number,
                           step_info=step_info, state=self.sample_set)

            # MCStep, i.e. actual sample move
            movepath = self._mover.move(self.sample_set, step=self.step)
            samples = movepath.results
            new_sampleset = self.sample_set.apply_samples(samples)
            elapsed_step = time.time() - time_start

            # TODO: we can save this with the MC steps for timing? The bit
            # below works, but is only a temporary hack
            setattr(movepath.details, "timing", elapsed_step)

            mcstep = MCStep(
                simulation=self,
                mccycle=self.step,
                previous=self.sample_set,
                active=new_sampleset,
                change=movepath
            )
            self._current_step = mcstep
            self.save_current_step()
            self.sample_set = new_sampleset

            # run after_step hooks
            hook_state = self.run_hooks('after_step', sim=self,
                                        step_number=step_number,
                                        step_info=step_info,
                                        # TODO: do we want the new state here
                                        # i.e. state after step as is right now
                                        # or do we want something similar to
                                        # shoot_snapshots where
                                        # old state == shooting snapshot
                                        state=self.sample_set,
                                        results=mcstep,
                                        hook_state=hook_state
                                        )
        # after simulation hooks
        self.run_hooks('after_simulation', sim=self)
