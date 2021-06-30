from __future__ import absolute_import
from nose.tools import (assert_equal, assert_not_equal, assert_almost_equal,
                        raises)
from nose.plugins.skip import Skip, SkipTest

import io
import sys
import time
import re
if ((sys.version_info.major < 3)
        or (sys.version_info.major == 3 and sys.version_info.minor < 4)):
    # re.fullmatch was introduced in py v3.4
    # but match should also do the trick
    re.fullmatch = re.match

import pytest
try:
    from unittest.mock import MagicMock
except ImportError:
    from mock import MagicMock

from ..test_helpers import (
    true_func, assert_equal_array_array, make_1d_traj, data_filename,
    assert_items_equal
)

import openpathsampling as paths

from openpathsampling.beta.hooks import *

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


class TestStorageHook(object):
    def setup(self):
        self.storage = MagicMock()
        self.empty_hook = StorageHook()
        self.hook = StorageHook(storage=self.storage,
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

    @pytest.mark.parametrize('step_num', [0, 5, 10])
    def test_after_step(self, step_num):
        self.hook.after_step(self.simulation, step_num, ('step', 'info'),
                             ('state'), "results", "hook_state")
        self.storage.save.assert_called_once()
        if step_num in [0, 10]:
            self.storage.sync_all.assert_called_once()
            # also check that we performed the sample sanity-ckeck
            self.simulation.sample_set.sanity_check.assert_called_once()

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
    # NOTE: we use regex to match slightly varying runtimes also
    # NOTE 2: waiting for 7 seconds is quite long, but the only way to get the
    # progress string for T > 1 min?
    def test_before_step(self, nn, elapsed):
        nn = nn  # MCstep counter during this simulation (reset to 0 every run)
        n_steps = 10  # Total number of MCSteps to do this simulation run
        # elapsed is total time elapsed since starting this simulation run
        step_info = nn, n_steps
        self.hook.before_simulation(self.simulation)  # start the hooks timer
        time.sleep(elapsed)
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
            assert re.fullmatch(
                    pattern=("Working on Monte Carlo cycle number 2\n"
                             + "Running for 0 seconds -  0.1[0-1] seconds per step\n"
                             + "Estimated time remaining: 0 seconds\n"),
                    string=contents)
        elif nn == 1 and elapsed == 1.:
            # should get progress string,
            # 10 s per step, 9 seconds (=9 steps) remaining
            assert re.fullmatch(
                    pattern=("Working on Monte Carlo cycle number 2\n"
                             + "Running for 1 second -  1.0[0-1] seconds per step\n"
                             + "Estimated time remaining: 9 seconds\n"),
                    string=contents)
        elif nn == 1 and elapsed == 7.:
            # should get progress string,
            # 7 s per step, 1:03 min (= 63 seconds (=9 steps)) remaining
            assert re.fullmatch(
                    pattern=("Working on Monte Carlo cycle number 2\n"
                             + "Running for 7 seconds -  7.0[0-1] seconds per step\n"
                             + "Estimated time remaining: 1 minute 3 seconds\n"),
                    string=contents)

    def test_after_simulation(self):
        # set step_number in PathSampling Mock
        self.simulation.configure_mock(step=2)
        self.hook.after_simulation(self.simulation, hook_state={})
        contents = self.stream.getvalue()
        assert contents == "DONE! Completed 2 Monte Carlo cycles.\n"
