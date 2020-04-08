import collections
import random
import numpy as np
from openpathsampling.volume import Volume
from openpathsampling.pathsimulators import PathSimulator
import openpathsampling as paths
from openpathsampling.netcdfplus import StorableNamedObject

import logging
logger = logging.getLogger(__name__)

AMSInfo = collections.namedtuple(
    'AMSInfo',
    ['initial_state', 'final_state', 'parametrized_volume']
)

class ParametrizedVolume(Volume):
    """Volume that can change by resetting some parameters.
    """
    def __init__(self, cv=None, cv_max=None):
        super(ParametrizedVolume, self).__init__()
        self._volume = None
        self.cv = cv
        self.cv_max = cv_max

    def set_parameters(self, parameters):
        """Use the given parameters to set the internal volume.

        Parameters
        ----------
        parameters : dict
            dictionary of parameters; depends on specific subclass

        See also
        --------
        unset_parameters
        """
        self._volume = self.volume_for_parameters(parameters)

    def unset_parameters(self):
        """Clears the previously used parameters.

        See also
        --------
        set_parameters
        """
        self._volume = None

    def volume_for_parameters(self, parameters):
        """Returns the volume associated with the given parameters.

        This is usually the primary method a subclass must implement.

        Parameters
        ----------
        parameters : dict
            dictionary of parameters; depends on specific subclass

        Returns
        -------
        :class:`.Volume` :
            volume associated with the given parameters
        """
        raise NotImplementedError()

    def __call__(self, snapshot):
        if self._volume:
            return self._volume(snapshot)
        else:
            raise RuntimeError("Cannot call parametrized volume before "
                               "setting parameters")

    def __str__(self):  # pragma: no cover
        return repr(self)  # str used in equality; default is wrong


class InterfaceSetParametrizedVolume(ParametrizedVolume):
    """Parameterized volume based on :class:`.VolumeInterfaceSet`

    Parameters
    ----------
    interface_set :class:`.VolumeInterfaceSet`
        interface set that can create volumes for this parameterized volume
    """
    def __init__(self, interface_set):
        super(InterfaceSetParametrizedVolume, self).__init__(
            cv=interface_set.cv,
            cv_max=interface_set.cv_max
        )
        self.interface_set = interface_set

    @classmethod
    def from_increasing_cv(cls, cv, min_val=float('-inf')):
        # this is hacky, but required because we infer direction of volume
        # interface set from number of volumes given (something to fix!)
        interface_set = paths.VolumeInterfaceSet(
            cv=cv,
            minvals=min_val,
            maxvals=[float('inf'), float('inf')]
        )
        return cls(interface_set)

    def volume_for_parameters(self, parameters):
        return self.interface_set.new_interface(**parameters)


class AMSInitialization(StorableNamedObject):
    """Intitial trajectories for adaptive multilevel splitting.
    """
    def __init__(self, initial_trajectories):
        self._initial_trajectories = initial_trajectories

    def to_dict(self):
        # default from_dict should work fine
        return {'initial_trajectories': self._initial_trajectories}

    @property
    def _default_initial_trajectories(self):
        return self._initial_trajectories

    @property
    def initial_trajectories(self):
        if self._initial_trajectories is None:
            self._initial_trajectories = self._default_initial_trajectories

        return self._initial_trajectories

class DynamicsAMSInitialization(AMSInitialization):
    """
    Use MD to run a long trajectory to turn into AMS initial trajectories.

    Parameters
    ----------
    initial_state : :class:`.Volume`
    final_state : :class:`.Volume`
    """
    class InitialTrajectoriesEnsemble(paths.WrappedEnsemble):
        """Ensemble to (fwd) generate initial AMS trajetories.
        """
        # TODO: we can add output to the can_append of this, like the
        # VisitAllStatesEnsemble
        def __init__(self, initial_state, final_state, n_trajectories):
            out_X = paths.AllOutXEnsemble(initial_state)
            in_X = paths.AllInXEnsemble(initial_state)
            sequence = ([paths.OptionalEnsemble(out_X)]
                        + [in_X, out_X] * n_trajectories
                        + [in_X & paths.LengthEnsemble(1)])
            ensemble = paths.SequentialEnsemble(sequence) \
                    & paths.AllOutXEnsemble(final_state)

            _cls = DynamicsAMSInitialization.InitialTrajectoriesEnsemble
            super(_cls, self).__init__(ensemble)

    def __init__(self, initial_state, final_state, initial_snapshot,
                 engine, n_trajectories=None, run=False):
        super(DynamicsAMSInitialization, self).__init__(
            initial_trajectories=None
        )
        self.initial_state = initial_state
        self.final_state = final_state
        self.initial_snapshot = initial_snapshot
        self.n_trajectories = n_trajectories
        self.engine = engine
        if run:
            self._initial_trajectories = self._default_initial_trajectories

    def generate_trajectories(self, initial_snapshot=None,
                              n_trajectories=None):
        """Run MD to actually create the trajectories
        """
        if n_trajectories is None:
            n_trajectories = self.n_trajectories

        if n_trajectories is None:
            raise ValueError("DynamicsAMSInitialization requires "
                             "n_trajectories, either in initialization or "
                             "when run() is called.")

        trajectories = []
        while len(trajectories) < n_trajectories:
            ensemble = self.InitialTrajectoriesEnsemble(
                initial_state=self.initial_state,
                final_state=self.final_state,
                n_trajectories=n_trajectories - len(trajectories)
            )
            traj = self.engine.generate(initial_snapshot,
                                        [ensemble.can_append])
            split_ens = paths.TISEnsemble(self.initial_state,
                                          self.final_state,
                                          self.initial_state)
            trajectories += split_ens.split(traj)
        assert len(trajectories) == n_trajectories, \
                "%d != %d" % (len(trajectories), n_trajectories)
        return trajectories

    @property
    def _default_initial_trajectories(self):
        return self.generate_trajectories(self.initial_snapshot,
                                          self.n_trajectories)


class AdaptiveMultilevelSplittingStepper(paths.SubPathMover):
    def __init__(self, initial_state, final_state, parametrized_volume,
                 engine):
        self.parametrized_volume = parametrized_volume
        self.ensemble = paths.TISEnsemble(initial_state, final_state,
                                          parametrized_volume)
        mover = paths.ForwardShootMover(
            ensemble=self.ensemble,
            selector=paths.InterfaceConstrainedSelector(parametrized_volume),
            engine=engine
        )
        super(AdaptiveMultilevelSplittingStepper, self).__init__(mover)

    def to_dict(self):
        dct = super(AdaptiveMultilevelSplittingStepper, self).to_dict()
        dct.update({
            'parametrized_volume': self.parametrized_volume,
            'ensemble': self.ensemble
        })
        return dct

    @classmethod
    def from_dict(cls, dct):
        obj = cls.__new__(cls)
        obj.mover = dct['mover']
        obj.parametrized_volume = dct['parametrized_volume']
        obj.ensemble = dct['ensemble']
        return obj

    @classmethod
    def from_AMS(cls, simulator):
        return cls(initial_state=simulator.initial_state,
                   final_state=simulator.final_state,
                   parametrized_volume=simulator.parametrized_volume,
                   engine=simulator.engine)

    def _select_min_lambda_trajectory(self, sample_set):
        # TODO: check with AMS guys to be sure that this is the correct way
        lambda_sample = collections.defaultdict(list)
        for samp in sample_set:
            lambda_val = self.parametrized_volume.cv_max(samp.trajectory)
            lambda_sample[lambda_val].append(samp)

        min_lambda = min(lambda_sample.keys())
        choice = random.choice(lambda_sample[min_lambda])
        return min_lambda, choice

    def move(self, sample_set):
        # select the trajectory to be changed
        (min_lambda, selected_sample) = \
            self._select_min_lambda_trajectory(sample_set)
        # select minimum max-lambda
        self.parametrized_volume.set_parameters({'lambda_i': min_lambda})
        # select random other traj
        all_replicas = set(sample_set.replica_list())
        replicas = list(all_replicas - {selected_sample.replica})
        if len(replicas) > 0:
            choice = int(np.random.choice(replicas))
        else:
            # special case for only 1 replica in sample set
            choice = int(selected_sample.replica)
        template_traj = sample_set[choice]
        # make sample set with other traj and selected repid
        new_starting_set = paths.SampleSet([paths.Sample(
            replica=selected_sample.replica,
            ensemble=selected_sample.ensemble,
            trajectory=template_traj.trajectory
        )])

        logger.info("Replica " + str(selected_sample.replica) + " at "
                    + str(min_lambda) + " replaced by shot from trajectory "
                    + "in replica " + str(choice))

        # run shooting move
        subchange = self.mover.move(new_starting_set)
        details = paths.Details(
            min_lambda=min_lambda,
            min_max_cv_sample=selected_sample,
            new_starting_set=new_starting_set
        )
        change = paths.SubMoveChange(
            subchange=subchange,
            mover=self,
            details=details
        )
        return change


class AdaptiveMultilevelSplitting(PathSimulator):
    """Main simulation class for adaptive multilevel splitting.
    """
    def __init__(self, storage, initial_state, final_state,
                 parametrized_volume, engine, initialization):
        super(AdaptiveMultilevelSplitting, self).__init__(storage)
        self.initial_state = initial_state
        self.final_state = final_state
        self.parametrized_volume = parametrized_volume
        self.engine = engine
        self.stepper = AdaptiveMultilevelSplittingStepper.from_AMS(self)
        self.initialization = initialization

        self.status_update_frequency = 1

    def _complete_count(self, sample_set):
        in_final_state = [self.final_state(sample.trajectory[-1])
                          for sample in sample_set]
        return sum(in_final_state)

    def run(self, complete_fraction=0.5, n_steps_max=None):
        initial_trajectories = self.initialization.initial_trajectories

        n_trans_stop = complete_fraction * len(initial_trajectories)
        self.sample_set = paths.SampleSet([
            paths.Sample(replica=repid,
                         ensemble=self.stepper.ensemble,
                         trajectory=traj)
            for (repid, traj) in enumerate(initial_trajectories)
        ])

        step = 0

        if self.storage:
            self.storage.steps.save(paths.MCStep(
                simulation=self,
                mccycle=step,
                active=self.sample_set,
                change=paths.AcceptedSampleMoveChange(self.sample_set.samples)
            ))

        continue_condition = lambda step, sample_set: (
            step < n_steps_max
            and self._complete_count(sample_set) < n_trans_stop
        )
        while continue_condition(step, self.sample_set):
            if step == 0 or (step + 1)  % self.status_update_frequency == 0:
                self.output(step, n_steps_max, complete_fraction,
                            self.sample_set, done=False)
            change = self.stepper.move(self.sample_set)
            self.sample_set = self.sample_set.apply_samples(change)
            mcstep = paths.MCStep(
                simulation=self,
                mccycle=step,
                active=self.sample_set,
                change=change
            )
            if self.storage:
                self.storage.steps.save(mcstep)
            step += 1

        self.output(step, n_steps_max, complete_fraction, self.sample_set,
                    done=True)

    def output(self, step, n_steps_max, complete_fraction, active, done):
        if done:
            out_str = "DONE: Ran {step} steps (maximum {n_steps_max})\n"
        else:
            out_str = "Running step {step} (maximum {n_steps_max})\n"
        out_str += "Max-max-lambda: {max_max}\n"
        out_str += "Min-max-lambda: {min_max}\n"
        lambda_vals = sorted([
            self.parametrized_volume.cv_max(samp.trajectory)
            for samp in active
        ])
        min_max = lambda_vals[0]
        max_max = lambda_vals[-1]
        string = out_str.format(step=step+1, n_steps_max=n_steps_max,
                                max_max=max_max, min_max=min_max)

        refresh = self.allow_refresh
        paths.tools.refresh_output(string,
                                   refresh=refresh,
                                   output_stream=self.output_stream)

    def analysis(self, steps=None):
        return AdaptiveMultilevelSplittingAnalysis.from_AMS(self, steps)


def DefaultAMS(storage, initial_state, final_state, cv, engine,
               initial_snapshot, n_initial_trajectories):
    initialization = DynamicsAMSInitialization(
        initial_state=initial_state,
        final_state=final_state,
        engine=engine,
        n_trajectories=n_initial_trajectories,
        initial_snapshot=initial_snapshot
    )
    parametrized_volume = \
        InterfaceSetParametrizedVolume.from_increasing_cv(cv)
    simulation = AdaptiveMultilevelSplitting(
        storage=storage,
        initial_state=initial_state,
        final_state=final_state,
        parametrized_volume=parametrized_volume,
        engine=engine,
        initialization=initialization
    )
    return simulation


class AdaptiveMultilevelSplittingAnalysis(StorableNamedObject):
    def __init__(self, initial_state, final_state, order_parameter, steps=None):
        self.initial_state = initial_state
        self.final_state = final_state
        self.order_parameter = order_parameter

    @classmethod
    def from_AMS(cls, simulation, steps=None):
        return cls(initial_state=simulation.initial_state,
                   final_state=simulation.final_state,
                   order_parameter=simulation.order_parameter,
                   steps=steps)

    def calculate_steps(self, steps):
        pass
