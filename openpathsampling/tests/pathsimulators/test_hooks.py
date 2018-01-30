from __future__ import absolute_import
from nose.tools import (assert_equal, assert_not_equal, assert_almost_equal,
                        raises)
from nose.plugins.skip import Skip, SkipTest

from ..test_helpers import (
    true_func, assert_equal_array_array, make_1d_traj, data_filename,
    assert_items_equal
)

import openpathsampling as paths

from openpathsampling.pathsimulators.hooks import *

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
        for step in range(n_steps):
            self.run_hooks('before_step', sim=self, step_number=step,
                           step_info=(step, n_steps), state=state)
            result = str(step)
            state += result
            self.run_hooks('after_step', sim=self, step_number=step,
                           step_info=(step, n_steps), state=state,
                           results=result, hook_state=hook_state)
        self.run_hooks('after_simulation', sim=self)

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

    def after_simulation(self, sim):
        self.finished_simulation = True

class StupiderHook(PathSimulatorHook):
    def __init__(self):
        self.foo = None

    def hook_method(self, sim, extra=None):
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
        extra = lambda sim: stupider.hook_method(sim, "foo")
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
    def _delete_files(self):
        pass

    def setup(self):
        pass

    def test_before_simulation(self):
        raise SkipTest

    def test_after_step(self):
        raise SkipTest

    def test_after_simulation(self):
        raise SkipTest

class TestShootFromSnapshotsOutputHook(object):
    def setup(self):
        pass

    def test_before_simulation(self):
        raise SkipTest

    def test_before_step(self):
        raise SkipTest
