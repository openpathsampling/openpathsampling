from __future__ import absolute_import
from nose.tools import (assert_equal, assert_not_equal, assert_almost_equal,
                        raises)
from nose.plugins.skip import Skip, SkipTest

import io
import sys
"""
import time
import re
if ((sys.version_info.major < 3)
        or (sys.version_info.major == 3 and sys.version_info.minor < 4)):
    # re.fullmatch was introduced in py v3.4
    # but match should also do the trick
    re.fullmatch = re.match
"""

import pytest
try:
    from unittest.mock import MagicMock, patch
except ImportError:
    from mock import MagicMock, patch

from ..test_helpers import (
    true_func, assert_equal_array_array, make_1d_traj, data_filename,
    assert_items_equal
)

import openpathsampling as paths

from openpathsampling.beta.hooks import *
from openpathsampling.beta.hooks import _self_or_sim_property_or_err

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)


class TrivialPathSimulator(paths.PathSimulator):
    """
    use this to check that the hooks get attached and run correctly
    """
    def run(self, n_steps):
        self.run_hooks('before_simulation', sim=self)
        state = ""
        hook_state = {}
        for step in range(n_steps):
            self.run_hooks('before_step', sim=self, step_number=step,
                           step_info=(step, n_steps), state=state)
            result = str(step)
            state += result
            hook_state = self.run_hooks('after_step', sim=self,
                                        step_number=step,
                                        step_info=(step, n_steps),
                                        state=state, results=result,
                                        hook_state=hook_state)
        self.run_hooks('after_simulation', sim=self, hook_state=hook_state)


class NoStepHookPathSimulator(paths.PathSimulator):
    """use this to check that only correct hooks get assigned"""
    hook_names = ['before_simulation', 'extra_hook']
    def run(self, n_steps):
        self.run_hooks('before_simulation', sim=self)
        self.run_hooks('extra_hook', sim=self, extra="foo")


class StupidHook(PathSimulatorHook):
    def __init__(self):
        self.began_simulation = False
        self.began_steps = 0
        self.finished_steps = 0
        self.finished_simulation = False

    def before_simulation(self, sim):
        self.began_simulation = True

    def before_step(self, sim, step_number, step_info, state):
        self.began_steps +=1

    def after_step(self, sim, step_number, step_info, state, results,
                   hook_state):
        self.finished_steps += 1

    def after_simulation(self, sim, hook_state):
        self.finished_simulation = True


class StupiderHook(PathSimulatorHook):
    def __init__(self):
        self.foo = None

    def hook_method(self, sim, hook_state=None, extra=None):
        if hook_state is None:
            hook_state = {}
        self.foo = extra


class TestPathSimulatorHook(object):
    # includes tests for attach_hook and run_hook
    def test_attach_hooks_class(self):
        hook = StupidHook()
        sim = TrivialPathSimulator(storage=None)
        sim.attach_hook(hook)
        for k in sim.hook_names:
            assert_equal(len(sim.hooks[k]), 1)

    def test_attach_hooks_method(self):
        hook_instance = StupiderHook()
        sim = NoStepHookPathSimulator(storage=None)
        sim.attach_hook(hook_instance.hook_method, hook_for='extra_hook')
        assert_not_equal(hook_instance.foo, "foo")
        sim.run(1)
        assert_equal(hook_instance.foo, "foo")

    def test_trivial_simulation_hooks(self):
        sim = TrivialPathSimulator(storage=None)
        stupid = StupidHook()
        stupider = StupiderHook()
        sim.attach_hook(stupid)
        extra = lambda sim, hook_state: stupider.hook_method(sim,
                                                             hook_state,
                                                             "foo")
        sim.attach_hook(extra, hook_for='after_simulation')
        # before running, check that state is as expected
        assert_equal(stupid.began_simulation, False)
        assert_equal(stupid.finished_simulation, False)
        assert_equal(stupid.began_steps, 0)
        assert_equal(stupid.finished_steps, 0)
        assert_not_equal(stupider.foo, "foo")
        # run
        sim.run(3)
        # check results
        assert_equal(stupid.began_simulation, True)
        assert_equal(stupid.finished_simulation, True)
        assert_equal(stupid.began_steps, 3)
        assert_equal(stupid.finished_steps, 3)
        assert_equal(stupider.foo, "foo")

    @raises(TypeError)
    def test_bad_extra_hook_method(self):
        sim = NoStepHookPathSimulator(storage=None)
        stupid = StupidHook()
        sim.attach_hook(stupid)

    @raises(TypeError)
    def test_bad_hook_name(self):
        sim = NoStepHookPathSimulator(storage=None)
        stupider = StupiderHook()
        sim.attach_hook(stupider.hook_method, "blah")

    def test_no_step_hook_simulation_hooks(self):
        sim = NoStepHookPathSimulator(storage=None)
        stupid = StupidHook()
        stupider = StupiderHook()
        sim.attach_hook(stupid.before_simulation,
                        hook_for='before_simulation')
        sim.attach_hook(stupider.hook_method, hook_for='extra_hook')
        assert_equal(len(sim.hooks), len(sim.hook_names))
        for hook in sim.hooks:
            assert_equal(len(sim.hooks[hook]), 1)

        assert_equal(stupid.began_simulation, False)
        assert_not_equal(stupider.foo, "foo")
        sim.run(1)
        assert_equal(stupid.began_simulation, True)
        assert_equal(stupider.foo, "foo")


class TestSelfOrSimProperty(object):
    def setup(self):
        self.attr_name = "test_attr"
        self.simulation = MagicMock(test_attr="sim_attr")
        self.sans_attr_sans_sim = MagicMock(_test_attr=None, _simulation=None)
        # Note: for the hooks test-attr is private, i.e. _test_attr
        self.with_attr_sans_sim = MagicMock(_test_attr="hook_attr")
        self.with_attr_with_sim = MagicMock(_test_attr="hook_attr",
                                            _simulation=self.simulation,
                                            )
        self.sans_attr_with_sim = MagicMock(_test_attr=None,
                                            _simulation=self.simulation,
                                            )
        # to test different names for the attributes on hook and sim
        self.attr_with_different_name_sans_sim = MagicMock(_attr="hook_attr")
        # set attr on hook to `None` to make sure we try to get it from sim
        self.attr_with_different_name_with_sim = MagicMock(
                                                    _attr=None,
                                                    _simulation=self.simulation,
                                                           )

    @pytest.mark.parametrize("hook_name", ["sans_attr_sans_sim",
                                           "with_attr_sans_sim",
                                           "with_attr_with_sim",
                                           "sans_attr_with_sim",
                                           "attr_with_different_name_with_sim",
                                           "attr_with_different_name_sans_sim",
                                           ])
    def test_self_or_sim_property_helper_func(self, hook_name):
        hook = getattr(self, hook_name)
        if hook_name == "sans_attr_sans_sim":
            with pytest.raises(SimulationNotFoundError,
                               match=("'" + self.attr_name + "' has not "
                                      + "been set and no hosting "
                                      + "simulation known to get "
                                      + "simulation." + self.attr_name + " ."),
                               ):
                _self_or_sim_property_or_err(hook, self.attr_name)
        elif hook_name == "with_attr_sans_sim":
            assert (_self_or_sim_property_or_err(hook, self.attr_name)
                    == "hook_attr")
        elif hook_name == "with_attr_with_sim":
            assert (_self_or_sim_property_or_err(hook, self.attr_name)
                    == "hook_attr")  # hook attribute takes precedence over sim
        elif hook_name == "sans_attr_with_sim":
            assert (_self_or_sim_property_or_err(hook, self.attr_name)
                    == "sim_attr")  # attr from sim, because hook-attr is not set
        # these two are the 'same' as the two above except for the differing attr names
        elif hook_name == "attr_with_different_name_sans_sim":
            assert (_self_or_sim_property_or_err(hook, "attr", self.attr_name)
                    == "hook_attr")
        elif hook_name == "attr_with_different_name_with_sim":
            assert (_self_or_sim_property_or_err(hook, "attr", self.attr_name)
                    == "sim_attr")


class TestStorageHook(object):
    def setup(self):
        self.storage = MagicMock(spec=["stash", "sync_all"])
        self.old_storage = MagicMock(spec=["save", "sync_all"])
        self.empty_hook = StorageHook()
        self.hook = StorageHook(storage=self.storage,
                                frequency=10)
        self.hook_old_storage = StorageHook(storage=self.old_storage,
                                            frequency=10)
        self.simulation = MagicMock(storage=self.storage,
                                    save_frequency=10)

    @pytest.mark.parametrize('hook', ['empty', 'std'])
    def test_before_simulation(self, hook):
        hook = {'empty': self.empty_hook,
                'std': self.hook}[hook]
        hook.before_simulation(self.simulation)
        assert hook.storage == self.storage
        assert hook.frequency == 10

    @pytest.mark.parametrize("hook_name", ["new_storage", "old_storage"])
    @pytest.mark.parametrize('step_num', [0, 5, 10])
    def test_after_step(self, step_num, hook_name):
        hook = {"new_storage": self.hook,
                "old_storage": self.hook_old_storage}[hook_name]
        storage = {"new_storage": self.storage,
                   "old_storage": self.old_storage}[hook_name]
        hook.after_step(self.simulation, step_num, ('step', 'info'),
                             ('state'), "results", "hook_state")
        if hook_name == "new_storage":
            storage.stash.assert_called_once()
        elif hook_name == "old_storage":
            storage.save.assert_called_once()
        if step_num in [0, 10]:
            storage.sync_all.assert_called_once()

    def test_after_simulation(self):
        self.hook.after_simulation(self.simulation, {})
        self.storage.sync_all.assert_called_once()


class TestShootFromSnapshotsOutputHook(object):
    def setup(self):
        if sys.version_info < (3,):
            self.stream = io.BytesIO()
        else:
            self.stream = io.StringIO()
        self.simulation = MagicMock(output_stream=self.stream,
                                    allow_refresh=False)
        self.empty_hook = ShootFromSnapshotsOutputHook()
        self.hook = ShootFromSnapshotsOutputHook(output_stream=self.stream,
                                                 allow_refresh=False)

    @pytest.mark.parametrize('hook_name', ['empty', 'std'])
    def test_before_simulation(self, hook_name):
        hook = {'empty': self.empty_hook,
                'std': self.hook}[hook_name]
        hook.before_simulation(self.simulation)
        assert hook.output_stream == self.stream
        assert hook.allow_refresh is False

    def test_before_step(self):
        step_info = (0, 100, 1, 10)  # snap_num, n_snap, shot_num, n_shots
        self.hook.before_step(self.simulation, step_number=2,
                              step_info=step_info, state=None)
        contents = self.stream.getvalue()
        assert contents == "Working on snapshot 1 / 100; shot 2 / 10"


class TestSampleSetSanityCheckHook(object):
    def setup(self):
        self.simulation = MagicMock(save_frequency=10)
        self.empty_hook = SampleSetSanityCheckHook()
        self.hook = SampleSetSanityCheckHook(frequency=10)

    @pytest.mark.parametrize('hook_name', ["empty", "std"])
    def test_before_simulation(self, hook_name):
        hook = {'empty': self.empty_hook,
                'std': self.hook}[hook_name]
        hook.before_simulation(self.simulation)
        assert hook.frequency == 10

    @pytest.mark.parametrize('step_num', [0, 5, 10])
    def test_after_step(self, step_num):
        self.hook.after_step(self.simulation, step_num, ('step', 'info'),
                             ('state'), "results", "hook_state")
        if step_num in [0, 10]:
            self.simulation.sample_set.sanity_check.assert_called_once()
        else:
            # check should only have been called for step_num % 10 == 0
            assert not self.simulation.sample_set.sanity_check.called


class TestLiveVisualizerHook(object):
    def setup(self):
        self.live_visualizer = MagicMock()
        self.simulation = MagicMock(live_visualizer=self.live_visualizer,
                                    status_update_frequency=2)
        self.empty_hook = LiveVisualizerHook()
        self.hook = LiveVisualizerHook(live_visualizer=self.live_visualizer,
                                       status_update_frequency=2)

    @pytest.mark.parametrize('hook_name', ['empty', 'std'])
    def test_before_simulation(self, hook_name):
        hook = {'empty': self.empty_hook,
                'std': self.hook}[hook_name]
        hook.before_simulation(self.simulation)
        assert hook.live_visualizer is self.live_visualizer
        assert hook.status_update_frequency == 2

    @pytest.mark.parametrize('step_num', [0, 1, 2])
    def test_after_step(self, step_num):
        self.hook.after_step(self.simulation, step_num, ('step', 'info'),
                             ('state'), "results", "hook_state")
        if step_num in [0, 2]:
            self.live_visualizer.draw_ipynb.assert_called_once()
        else:
            # should only visualize when step_num % 2 == 0
            assert not self.live_visualizer.draw_ipynb.called


class TestPathSamplingOutputHook(object):
    def setup(self):
        if sys.version_info < (3,):
            self.stream = io.BytesIO()
        else:
            self.stream = io.StringIO()
        self.simulation = MagicMock(output_stream=self.stream,
                                    allow_refresh=False)
        self.empty_hook = PathSamplingOutputHook()
        self.hook = PathSamplingOutputHook(output_stream=self.stream,
                                           allow_refresh=False,
                                           status_update_frequency=1)

    @pytest.mark.parametrize('hook_name', ['empty', 'std'])
    def test_before_simulation(self, hook_name):
        hook = {'empty': self.empty_hook,
                'std': self.hook}[hook_name]
        hook.before_simulation(self.simulation)
        assert hook.output_stream == self.stream
        assert hook.allow_refresh is False

    @pytest.mark.parametrize('nn, elapsed',
                             [(0, 0.1), (1, 0.1), (1, 1.), (1, 7.)]
                             )
    def test_before_step(self, nn, elapsed):
        nn = nn  # MCstep counter during this simulation (reset to 0 every run)
        n_steps = 10  # Total number of MCSteps to do this simulation run
        # elapsed is total time elapsed since starting this simulation run
        step_info = nn, n_steps
        with patch("openpathsampling.beta.hooks.time.time",
                   new=MagicMock(side_effect=[0, elapsed]),
                   ):
            self.hook.before_simulation(self.simulation)  # start the hooks timer
            self.hook.before_step(self.simulation, step_number=2,
                                  step_info=step_info, state=None)
        contents = self.stream.getvalue()
        if nn == 0:
            # we should get the start simulation string
            assert contents == ("Working on Monte Carlo cycle number 2\n"
                                + "Starting simulation...\n"
                                + "Working on first step\n")
        elif nn == 1 and elapsed == 0.1:
            # should get progress string,
            # 0.10 s per step, 0.9 seconds (=9 steps)) remaining
            assert (contents == "Working on Monte Carlo cycle number 2\n"
                                 + "Running for 0 seconds -  0.10 seconds per step\n"
                                 + "Estimated time remaining: 0 seconds\n"
                    )
        elif nn == 1 and elapsed == 1.:
            # should get progress string,
            # 1.00 s per step, 9 seconds (=9 steps) remaining
            assert (contents == "Working on Monte Carlo cycle number 2\n"
                                + "Running for 1 second -  1.00 seconds per step\n"
                                + "Estimated time remaining: 9 seconds\n"
                    )
        elif nn == 1 and elapsed == 7.:
            # should get progress string,
            # 7.00 s per step, 1:03 min (= 63 seconds (=9 steps)) remaining
            assert (contents == "Working on Monte Carlo cycle number 2\n"
                                + "Running for 7 seconds -  7.00 seconds per step\n"
                                + "Estimated time remaining: 1 minute 3 seconds\n"
                    )

    def test_after_simulation(self):
        # set step_number in PathSampling Mock
        self.simulation.configure_mock(step=2)
        self.hook.after_simulation(self.simulation, hook_state={})
        contents = self.stream.getvalue()
        assert contents == "DONE! Completed 2 Monte Carlo cycles.\n"


class TestGraciousKillHook(object):
    def setup(self):
        self.final_call = MagicMock()
        self.nocall_final_call = MagicMock()
        self.simulation = MagicMock(storage=MagicMock(close=MagicMock(),
                                                      sync_all=MagicMock()
                                                      )
                                    )
        self.nocall_simulation = MagicMock(storage=MagicMock(close=MagicMock(),
                                                             sync_all=MagicMock()
                                                             )
                                           )
        self.kill_hook = GraciousKillHook("0 hours 20 second",
                                          final_call=self.final_call
                                          )
        self.nokill_hook = GraciousKillHook("10 hours 2 seconds",
                                            final_call=self.nocall_final_call
                                            )

    def test_kill(self):
        s_num = 1  # step after which the simulation is killed
        with patch("openpathsampling.beta.hooks.time.time",
                   new=MagicMock(side_effect=[0, 10]),
                   ):
            # start the timer
            self.kill_hook.before_simulation(self.simulation)
            with pytest.raises(GraciousKillError):
                self.kill_hook.after_step(sim=self.simulation,
                                          step_number=s_num,
                                          step_info=(0, 200),
                                          state=None,
                                          results=None,
                                          hook_state=None
                                          )
        # now assert that everything got called
        self.simulation.storage.sync_all.assert_called_once()
        self.simulation.storage.close.assert_called_once()
        self.final_call.assert_called_once_with(s_num)

    def test_nokill(self):
        s_num = 1  # step after which the simulation is not killed
        with patch("openpathsampling.beta.hooks.time.time",
                   new=MagicMock(side_effect=[0, 10]),
                   ):
            # start the timer
            self.nokill_hook.before_simulation(self.simulation)
            self.nokill_hook.after_step(sim=self.simulation,
                                        step_number=s_num,
                                        step_info=(0, 200),
                                        state=None,
                                        results=None,
                                        hook_state=None
                                        )
        # now assert that everything got called
        self.nocall_simulation.storage.sync_all.assert_not_called()
        self.nocall_simulation.storage.close.assert_not_called()
        self.nocall_final_call.assert_not_called()
